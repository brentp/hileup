import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pysam
import cyvcf2
import chileup
import time

config = chileup.Config(tags=[], track_read_names=False,
            track_base_qualities=False, track_mapping_qualities=False,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)

bam = pysam.AlignmentFile(sys.argv[1], "rb",
                          reference_filename="/data/human/g1k_v37_decoy.fa")

vcf = cyvcf2.VCF(sys.argv[2])

t0 = time.time()

for v in vcf("1"):
    h = chileup.pileup(bam, v.CHROM, v.start, config)
    print(repr(v).strip(), h.bases)

print("hileup time: %.2f", time.time() - t0)

# now do pysam pileup.
