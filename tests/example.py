import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pysam
import cyvcf2
import chileup

config = chileup.Config(tags=[], track_read_names=True,
            track_base_qualities=True, track_mapping_qualities=True,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)

bam = pysam.AlignmentFile(sys.argv[1], "rb",
                          reference_filename="/data/human/g1k_v37_decoy.fa")

vcf = cyvcf2.VCF(sys.argv[2])

for v in vcf:
    h = chileup.pileup(bam, v.CHROM, v.start, config)
    print(repr(v).strip(), h.bases)
