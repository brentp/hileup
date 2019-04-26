import strutils
import tables

type base_strand* = object
  reverse_strand* {.bitsize: 1.}: uint8
  base* {.bitsize: 7.}: uint8

proc `$`(b:base_strand): string =
  if b.reverse_strand == 1:
    return $(b.base.char.toLowerAscii)
  return $(b.base.char)

type insertion* = object
  index*: uint32 ## index in the hile object to get qname
  len*: uint32
  sequence: string

type deletion* = object
  index*: uint32 ## index in the hile object to get qname
  len*: uint16

type HileSite* = object
  ## a hilesite holds all the info for a particular genomic position
  qname_idxs*: seq[uint32]
  base_strands*: seq[base_strand]
  bqs*: seq[uint8]
  mqs*: seq[uint8]

type QnameAt = object
  qname_i: uint32
  position: uint32

type Qnames = object
  qnames_u2s: TableRef[uint32, string]
  qnames_s2u: TableRef[string, uint32]

proc add(q:Qnames, qname:string): uint32 =
  if qname in q.qnames_s2u:
    return q.qnames_s2u[qname]
  result = q.qnames_s2u.len.uint32
  q.qnames_s2u[qname] = result
  q.qnames_u2s[result] = qname

template get(q:Qnames, qname: string): uint32 =
  q.qname_s2u[qname]

template get(q:Qnames, u: uint32): string =
  q.qname_u2s[u]

type Hiles* = object
  ## so hiles[i] is the info for genomic position i+refOffset
  ## qnames are optionally stored externally
  refOffset*: uint32
  # allow finding id by name and name by id. since every base in the read has a hilesite
  # we avoid storing the read-name more than needed.
  qnames: Qnames
  # a qnameat allows us to get the exact base, read combination where the
  # insertion /deletion occurred.
  insertions: TableRef[QnameAt, seq[insertion]]
  deletions: TableRef[QnameAt, seq[deletion]]
  hiles*: seq[HileSite]
