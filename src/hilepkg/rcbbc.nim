import argparse
import strformat
import strutils
import sets
import hts/bam
import hts/bgzf/bgzi
import lapper

include ./utils

type pair = object
    start: int
    stop: int

iterator gen_start_ends(c: Cigar, ipos: int): pair {.inline.} =
  # generate start, end pairs given a cigar string and a position offset.
  if c.len == 1 and c[0].op == CigarOp.match:
    yield pair(start:ipos, stop:ipos + c[0].len)
  else:
    var pos = ipos
    var last = pair(start:pos, stop: -1)
    var con: Consume
    for op in c:
      con = op.consumes
      if not con.reference:
        continue
      var olen = op.len
      if con.query:
        if pos != last.stop:
          if last.stop != -1:
            yield last
            last = pair(start: pos, stop: -1)
        last.stop = pos + olen
      pos += olen
    if last.stop != -1:
      yield last

proc output(bed:TableRef[string, Lapper[region]], order: seq[string], output_prefix: string) =
  var out_path = output_prefix & ".counts.bed.gz"
  var bgz = wopen_bgzi(out_path, 0, 1, 2, true)
  stderr.write_line &"[rcbbc] writing bgzipped long format file to {out_path}"
  for chrom in order:
    for reg in bed[chrom]:
      for bc, cnt in reg.barcodes:
        discard bgz.write_interval(&"{chrom}\t{reg.start}\t{reg.stop}\t{bc}\t{cnt}", chrom, reg.start, reg.stop)
  discard bgz.close()

proc output_sparse(bed:TableRef[string, Lapper[region]], order: seq[string], white_list: OrderedSet[string], output_prefix: string) =
  var out_path = output_prefix & ".mtx"
  var fh: File
  if not open(fh, out_path, fmWrite):
    quit "couldn't open sparse file:" & out_path

  var fhcells: File
  if not open(fhcells, out_path & ".cells", fmWrite):
    quit "couldn't open sparse cells label file:" & out_path & ".cells"

  var fhsites: File
  if not open(fhsites, out_path & ".sites", fmWrite):
    quit "couldn't open sparse cells sites file:" & out_path & ".sites"

  var sites_order = initTable[string, int](8192)
  var cells_order = initTable[string, int](8192)
  if white_list.len > 0:
    for bc in white_list:
      cells_order[bc] = cells_order.len
      fhcells.write_line bc

  var nsites = 0
  var nvalues = 0
  for chrom in order:
    nsites += bed[chrom].len
    for reg in bed[chrom]:
      if not sites_order.hasKeyOrPut(reg.id, sites_order.len):
        fhsites.write_line reg.id
      nvalues.inc(reg.barcodes.len)
      for bc, cnt in reg.barcodes:
        if white_list.len == 0 and not cells_order.hasKeyOrPut(bc, cells_order.len):
          fhcells.write_line bc

  fh.write("%%MatrixMarket matrix coordinate integer general\n")
  fh.write(&"{nsites} {cells_order.len} {nvalues}\n")
  for chrom in order:
    for reg in bed[chrom]:
      for bc, cnt in reg.barcodes:
        fh.write(&"{sites_order[reg.id]+1} {cells_order[bc]+1} {cnt}\n")

  fh.close()
  stderr.write_line &"[rcbbc] wrote a {sites_order.len}x{cells_order.len} sparse matrix to {out_path} with {out_path}.sites and {out_path}.cells for row and column labels"
  stderr.write_line &"""[rcbbc] read this in R with: library(Matrix); m = readMM("{out_path}"); rownames(m)=readLines("{out_path}.sites"); colnames(m)=readLines("{out_path}.cells")"""
  stderr.write_line &"""[rcbbc] read this in python with: from scipy import io as sio; import pandas as pd; A = sio.mmread("{out_path}"); df = pd.DataFrame(A.toarray(), index=[l.strip() for l in open("{out_path}.sites")], columns=[l.strip() for l in open("{out_path}.cells")])"""
  fhsites.close()
  fhcells.close()


proc output_qc(bed: TableRef[string, Lapper[region]], output_prefix: string, counts_by_bc_total: CountTable[string], counts_by_bc_target: CountTable[string]) =
  # columns are:
  # 1. cell-barcode
  # 2. total reads for this bar code.
  # 3. reads from this barcode within the BED regions.
  # 4. proportion of reads for this barcode  within the bed regions (FRiP)
  var out_path = output_prefix & ".cb-counts.txt"
  stderr.write_line &"[rcbbc] writing reads-per-cell and FRiP to {out_path}"
  var fh:File
  if not open(fh, out_path, fmWrite):
    quit "couldn't open path for cell-barcodes at:" & out_path

  fh.write("#barcode\ttotal_reads\ton_target_reads\tfraction_on_target\n")
  for bc in counts_by_bc_total.keys:
    var tot = counts_by_bc_total[bc]
    var ont = counts_by_bc_target.getOrDefault(bc)
    fh.write(&"{bc}\t{tot}\t{ont}\t{ont.float / tot.float:.3f}\n")
  fh.close()

proc rcbbc(ibam: Bam, bed: TableRef[string, Lapper[region]], white_list: OrderedSet[string], tagid: string, exclude_flags: uint16, min_mapping_quality: uint8, output_prefix:string) =

  var last_tid = -1
  var last_tree: Lapper[region]
  var order: seq[string]
  doAssert tagid.len == 2
  var total = 0
  var notag = 0
  var added = 0
  var non_white_list = 0
  var hasWhiteList = white_list.len > 0

  var counts_by_bc_total = initCountTable[string](8192)
  var counts_by_bc_target = initCountTable[string](8192)
  var reported = false


  for aln in ibam:
    total.inc
    if aln.mapping_quality < min_mapping_quality: continue
    if (aln.flag and exclude_flags) != 0: continue

    if aln.tid != last_tid:
      try:
        last_tree = bed[aln.chrom]
        order.add(aln.chrom)
      except KeyError:
        if not reported:
          stderr.write_line &"[rcbbc] chromosome {aln.chrom} not found in intervals in bed file (this is the last warning related to missing chroms)"
          reported = true
        last_tid = aln.tid
        var x:seq[region]
        last_tree = lapify(x)

      last_tid = aln.tid

    var bc = tag[cstring](aln, tagid)
    var btag: string
    if bc.isNone:
      btag = "notag"
      notag += 1
    else:
      btag = $bc.get
    if hasWhiteList and btag notin white_list:
      non_white_list.inc
      continue

    var on_target = false
    counts_by_bc_total.inc(btag)
    for p in gen_start_ends(aln.cigar, aln.start):
      last_tree.each_seek(p.start, p.stop, proc(r:region) = (r.barcodes.inc(btag); on_target = true))
    if on_target:
      counts_by_bc_target.inc(btag)
    added.inc

  stderr.write_line(&"[rcbbc] added {added} alignments out of {total}; of those {notag} had no tag for {tagid}")
  if hasWhiteList:
    stderr.write_line(&"[rcbbc] {non_white_list} alignments skipped because they were not present in the specified white-list")

  bed.output(order, output_prefix)
  bed.output_qc(output_prefix, counts_by_bc_total, counts_by_bc_target)
  bed.output_sparse(order, white_list, output_prefix)

proc read_white_list(path: string): OrderedSet[string] =
  result = initOrderedSet[string](8192)
  if path == "": return
  for l in path.lines:
    result.incl(l.strip())

proc main() =

  var p = newParser("rcbbc"):
    help("read-count by barcode")
    option("-f", "--fasta", help="path to fasta. required only for CRAM")
    option("-e", "--exclude", help="exclude alignments with any of these bits set", default="1796")
    option("-m", "--min-mapping-quality", help="exclude alignments with a mapping-quality below this", default="1")
    option("-w", "--white-list", help="line-delimited list of barcodes to include. default is to include any seen in the bam")
    option("-t", "--tag", default="CB", help="tag on which to partition read-counts")
    option("-p", "--prefix", default="rcbbc", help="output prefix for files written by rcbbc")
    arg("bed")
    arg("bam")

  try:
    if p.parse().help:
      quit 0
  except UsageError:
    echo p.help
    quit 2

  var opts = p.parse()
  if opts.help:
    quit 0
  var fai_path : cstring = nil
  if opts.fasta != "":
    fai_path = opts.fasta.cstring

  var bed = read_bed(opts.bed)
  var ibam: Bam
  if not open(ibam, opts.bam, fai=fai_path, threads=2):
    quit "[rcbbc] couldn't open bam file:" & opts.bam
  if opts.tag.len != 2:
    echo p.help()
    quit "[rcbbc] tag must be length 2"

  var white_list = read_white_list(opts.white_list)
  if opts.prefix == "rcbbc":
    stderr.write_line &"[rcbbc] outputting files with prefix {opts.prefix}; use --prefix to change this"

  ibam.rcbbc(bed, white_list, opts.tag, parseInt(opts.exclude).uint16, parseInt(opts.min_mapping_quality).uint8, opts.prefix)


when isMainModule:
  main()
