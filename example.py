import pysam
import chileup

bam = pysam.AlignmentFile("tests/three.bam", "rb")

config = chileup.Config(tags=[], track_read_names=True,
        track_base_qualities=True, track_mapping_qualities=True,
        exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
        min_base_quality=10, min_mapping_quality=10)


h = chileup.pileup(bam, "1", 1585270, config)

print(h.bases)
print(h.read_names)
print(h.bqs)
print(h.mqs)
print(h.deletions)
print(h.insertions)
print(h.tags)
