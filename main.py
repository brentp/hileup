import pysam
import chileup

bam = pysam.AlignmentFile("tests/three.bam", "rb")

# setting track_xxx to False will speed the hileup as less copying and and data access is required.
config = chileup.Config(tags=[], track_read_names=True,
        track_reads=True,
        track_base_qualities=True, track_mapping_qualities=True,
        exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
        min_base_quality=10, min_mapping_quality=10)


# We can ignore all reads with 'C' at this base (for example to get variant-only reads)
h = chileup.pileup(bam, "1", 1585270, config, 'C')

print(h.bases) # 'tt'
print(h.read_names) # [b'A00227:74:HCWC7DSXX:1:1269:13449:13855', b'A00227:74:HCWC7DSXX:1:2426:7157:15483']
print(h.bqs) # [37, 37]  numpy array that is a view into underlying data.
print(h.mqs) # [60, 60] numpy view.
print(h.deletions) # [(0, 8) (1, 8)] numpy view

# NOTE: if you're needing this, it might be simpler to use pysam pileup.
reads = h.reads(bam.header) # [<pysam.libcalignedsegment.AlignedSegment object at 0x7f3cd652bca0>, <pysam.libcalignedsegment.AlignedSegment object at 0x7f3cd506e500>]
read_positions = [chileup.query_pos(r) for r in reads]
print(read_positions) # [69, 69]]
print("query sequence:", [read.query_sequence[p] for (read, p) in zip(reads, read_positions)]) # 'TT' matches bases above.





# the insertions and deletions have a `.index` property that can be used
# to access the read-names, tags, etc that are associated with the indel event.
for ins in h.insertions: # copy of the data.
   print(h.read_names[ins.index], h.tags[ins.index], ins.sequence, ins.len)

print('tags:', h.tags) # copy.
