import argparse
import strformat
import hts/bam
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

proc output(bed:TableRef[string, Lapper[region]], order: seq[string]) =
  for chrom in order:
    for reg in bed[chrom]:
      for bc, cnt in reg.barcodes:
        echo &"{chrom}\t{reg.start}\t{reg.stop}\t{bc}\t{cnt}"

proc rcbbc(ibam: Bam, bed: TableRef[string, Lapper[region]], exclude_flags: uint16, min_mapping_quality: uint8) =

  var last_tid = -1
  var last_tree: Lapper[region]
  var order: seq[string]

  for aln in ibam:
    if aln.mapping_quality < min_mapping_quality: continue
    if (aln.flag and exclude_flags) != 0: continue

    if aln.tid != last_tid:
      try:
        last_tree = bed[aln.chrom]
        order.add(aln.chrom)
      except KeyError:
        stderr.write_line &"[rcbbc] chromosome {aln.chrom} not found in intervals in bed file"
        last_tid = aln.tid
        var x:seq[region]
        last_tree = lapify(x)

      last_tid = aln.tid

    var bc = tag[cstring](aln, "CB")
    if bc.isNone:
      continue

    for p in gen_start_ends(aln.cigar, aln.start):
      last_tree.each_seek(p.start, p.stop, proc(r:region) = r.barcodes.inc($bc.get))

  bed.output(order)


proc main() =

  var p = newParser("rcbbc"):
    option("-f", "--fasta", help="path to fasta. required only for CRAM")
    option("-e", "--exclude", help="exclude alignments with any of these bits set", default="1796")
    option("-m", "--min-mapping-quality", help="exclude alignments with a mapping-quality below this", default="1")
    option("-t", "--tag", default="CB", help="tag on which to partigion read-counts")
    flag("-w", "--wide", help="output counts in wide form (1 column per barcode) default is long form")
    arg("bed")
    arg("bam")

  var opts = p.parse()
  var fai_path : cstring = nil
  if opts.fasta != "":
    fai_path = opts.fasta.cstring


  var bed = read_bed(opts.bed)
  var ibam: Bam
  if not open(ibam, opts.bam, fai=fai_path, threads=2):
    quit "[rcbbc] couldn't open bam file:" & opts.bam


  ibam.rcbbc(bed, parseInt(opts.exclude).uint16, parseInt(opts.min_mapping_quality).uint8)


when isMainModule:
  main()
