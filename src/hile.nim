import hts
import strutils
import sets

type Config* = object
  MinMappingQuality*: uint8
  MinBaseQuality*: uint8
  ExcludeFlags*: uint16
  IncludeFlags*: uint16
  TrackReadNames*: bool
  TrackBaseQualities*: bool
  TrackMappingQualities*: bool
  ParitionTag*: seq[array[2, char]]

type insertion* = object
  index: uint16 ## index in the hile object to get qname
  sequence: string

type deletion* = object
  index*: uint16 ## index in the hile object to get qname
  len*: uint16

type basestrand* = object
  reverse_strand* {.bitsize: 1.}: uint8
  base* {.bitsize: 7.}: uint8

proc `$`(b:basestrand): string =
  if b.reverse_strand == 1:
    return $(b.base.char.toLowerAscii)
  return $(b.base.char)

type Hile* = object
  pos*: int32
  reference_base*: char
  bases*: seq[basestrand]
  bqs*: seq[uint8]
  mqs*: seq[uint8]
  read_names*: seq[string]
  ins*: seq[insertion]
  del*: seq[deletion]
  tag*: array[2, char]

proc size(h:Hile): int =
  result = sizeof(h) + sizeof(char) * h.bases.len + sizeof(uint8) * (h.mqs.len + h.bqs.len)
  for q in h.read_names:
    result += q.len
  result += sizeof(deletion) * h.del.len
  for ins in h.ins:
    result += sizeof(ins) + ins.sequence.len

proc primary(aln: Record): bool {.inline.} =
  return not (aln.flag.secondary or aln.flag.supplementary)

proc hileup*(bam:Bam, chrom: string, position:int, reference: Fai, cfg:Config, result: var Hile) {.inline.} =

  var overlap_lookup = initHashSet[string]()
  result.pos = position.int32
  if reference != nil:
    result.reference_base = reference.get(chrom, position, position + 1)[0]
  var s = ""
  for aln in bam.query(chrom, position, position + 1):
    if aln.mapping_quality < cfg.MinMappingQuality: continue
    if cfg.IncludeFlags != 0 and (cfg.IncludeFlags and aln.flag.uint16) != cfg.IncludeFlags: continue
    if cfg.ExcludeFlags != 0 and (cfg.ExcludeFlags and aln.flag.uint16) != 0: continue

    var
      r_off = aln.start
      q_off = 0
      skip_last = false # don't append INS if preceding base was skipped.

    for op in aln.cigar:
      if r_off == position + 1 and not skip_last:
        if op.op == CigarOp.deletion:
          result.del.add(deletion(index: result.bases.high.uint16, len: op.len.uint16))
        elif op.op == CigarOp.insert:
          discard aln.sequence(s)
          result.ins.add(insertion(index: result.bases.high.uint16, sequence: s[q_off..<q_off + op.len]))

      if r_off > position: break
      skip_last = false

      let cons = op.consumes
      if cons.query:
        q_off += op.len
      if cons.reference:
        r_off += op.len

      if r_off < position: continue

      if not (cons.query and cons.reference):
        continue

      var over = r_off - position
      if over > q_off: break
      if over < 0: continue

      if aln.tid == aln.mate_tid and overlap_lookup.len > 0 and not overlap_lookup.missingOrExcl(aln.qname):
        continue

      if cfg.MinBaseQuality > 0'u8 or cfg.TrackBaseQualities:
        var bq = aln.base_quality_at(int(q_off - over))
        if bq < cfg.MinBaseQuality:
          skip_last = true
          continue
        if cfg.TrackBaseQualities:
          result.bqs.add(bq)

      var c = aln.base_at(int(q_off - over))
      var bs = basestrand(base: c.uint8)
      if aln.flag.reverse:
        bs.reverse_strand = 1
      result.bases.add(bs)

      if cfg.TrackReadNames:
        result.read_names.add(aln.qname)
      if cfg.TrackMappingQualities:
        result.mqs.add(aln.mapping_quality)

      if not aln.primary or aln.stop <= aln.mate_pos or aln.mate_pos > position or aln.tid != aln.mate_tid or aln.start > aln.mate_pos:
        continue
      overlap_lookup.incl(aln.qname)

proc hileup*(bam:Bam, chrom: string, position:int, reference: Fai, cfg:Config): Hile =
    hileup(bam, chrom, position, reference, cfg, result)

when isMainModule:
  import hts/private/hts_concat
  #import nimprof

  import os
  var b:Bam
  var fai:Fai
  if not open(b, paramStr(1), index=true, fai = paramStr(2)):
    quit "bad bam"
  if not open(fai, paramStr(2)):
    quit "bad fai"

  var cfg = Config(MinBaseQuality: 9, MinMappingQuality: 9, ExcludeFlags: BAM_FUNMAP or BAM_FSECONDARY or BAM_FQCFAIL or BAM_FDUP)
  cfg.TrackBaseQualities = true
  cfg.TrackMappingQualities = true
  cfg.TrackReadNames = true
  #for i in 1584850..1586144:
  #for i in 1585268..1585280:
  # samtools mpileup -B -q 10 -Q 10 -r 1:20101-20101 --reference /data/human/g1k_v37_decoy.fa /data/human/hg002.cram
  for i in 20100..20100:
    echo i
    var h = b.hileup("1", i, fai, cfg)
    echo h
    if h.bases.len == 0: continue
    echo h.pos, " ", h.bases.len
