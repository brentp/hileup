import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pysam
import chileup

config = chileup.Config(tags=[], track_read_names=True,
            track_reads=True,
            track_base_qualities=True, track_mapping_qualities=True,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)

def test_off_by_one():
    bam = pysam.AlignmentFile("tests/soft.bam", "rb")
    pos = 10080
    h = chileup.pileup(bam, "1", pos, config)
    assert h.bases == "", h.bases

def test_insertion():
    bam = pysam.AlignmentFile("tests/three.bam", "rb")
    pos = 1585270
    h = chileup.pileup(bam, "1", pos, config)
    assert h.bases == 'tt', h.bases
    assert len(h.deletions) == 2

def test_reads():
    bam = pysam.AlignmentFile("tests/three.bam", "rb")
    pos = 1585270
    h = chileup.pileup(bam, "1", pos, config)
    assert h.bases == 'tt', h.bases
    assert len(h.deletions) == 2
    print(h.reads(bam.header))


def main(bam, config):

    actg = set("ATGCatgc")
    for pos in range(1585202, 1585202 + 150):
        h = chileup.pileup(bam, "1", pos, config)
        if len(h.insertions) > 0:
            print(pos, len(h.bases), h.insertions)
        #print(pos, h.bases, len(h.bqs), len(h.deletions), len(h.insertions))
        if set(h.bases) - actg:
            print(pos, h.bases, set(h.bases) - actg)
            break
        h.rbp(h)

if __name__ == "__main__":
    bam = pysam.AlignmentFile("tests/three.bam")

    config = chileup.Config(tags=[], track_read_names=True,
            track_base_qualities=True, track_mapping_qualities=True,
            track_reads=True,
            exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
            min_base_quality=10, min_mapping_quality=10)

    pos = 1585202
    h = chileup.pileup(bam, "1", pos, config)
    print(h)
    #print(pos, h.bases)
    test_off_by_one()
    print("off by one")
    test_insertion()
    print("insertion")

    test_reads()

    main(bam, config)
