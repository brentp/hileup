# distutils: sources = hile.c
# distutils: include_dirs = .

from pysam.libcalignmentfile cimport AlignmentFile, AlignmentHeader
from pysam.libcalignedsegment cimport AlignedSegment
from libc.stdlib cimport malloc
from collections import Counter
cimport numpy as np
import numpy as np
cdef dt = np.dtype([('index', np.uint32, 1), ('len', np.uint32, 1)])


cdef inline dict decode_insertion(dict ins):
    ins["sequence"] = ins["sequence"].decode()
    return ins

cdef class HileUp:
    cdef hile *c

    #@property
    #def reference_base(self):
    #    return chr(self.c.reference_base)

    @property
    def bases(HileUp self):
        cdef int i
        cdef char* s = <char *>malloc((self.c.n+1) * sizeof(char))
        for i in range(self.c.n):
            s[i] = self.c.bases[i].base + (self.c.bases[i].reverse_strand * 32)
        s[self.c.n] = 0
        pys = s[:self.c.n].decode()
        free(s)
        return pys

    def rbp(HileUp self, HileUp other):
        cdef dict self_lookup = {}
        cdef int i
        for i in range(self.c.n):
            self_lookup[self.c.read_names[i]] = chr(self.c.bases[i].base)
        result = {}
        for i in range(other.c.n):
            if other.c.read_names[i] not in self_lookup:
                continue
            key = self_lookup[other.c.read_names[i]] + chr(other.c.bases[i].base)
            if not key in result: result[key] = 0
            result[key] += 1
        return result

    @property
    def tags(self):
        cdef int i
        if self.c.tags == NULL: return []
        return [self.c.tags[i] for i in range(self.c.n)]

    @property
    def bqs(self):
        if self.c.bqs == NULL: return []
        return np.asarray(<np.uint8_t[:self.c.n]>self.c.bqs, dtype=np.uint8)

    @property
    def mqs(self):
        if self.c.mqs == NULL: return []
        return np.asarray(<np.uint8_t[:self.c.n]>self.c.mqs, dtype=np.uint8)

    @property
    def read_names(self):
        if self.c.read_names == NULL: return []

        cdef int i
        result = []
        for i in range(0, self.c.n):
            result.append(self.c.read_names[i][:])
        return result

    def reads(self, AlignmentHeader hdr):
        if self.c.reads == NULL: return []
        cdef int i
        result = []
        cdef AlignedSegment dest
        for i in range(0, self.c.n):
            dest = AlignedSegment.__new__(AlignedSegment)
            dest._delegate = bam_dup1(self.c.reads[i])
            dest.header = hdr
            result.append(dest)
        return result

    @property
    def deletions(self):
        """return a view of the data about deletions."""
        if self.c.n_deletions == 0: return []
        return np.asarray(<hile_deletion_t[:self.c.n_deletions]>self.c.deletions, dtype=dt)

    @property
    def insertions(self):
        "return the information about the insertions. This is copy of the underlying data."
        if self.c.n_insertions == 0: return []
        return [decode_insertion(self.c.insertions[i]) for i in range(self.c.n_insertions)]

    @property
    def size(self):
        return self.c.n

    def __dealloc__(self):
        hile_destroy(self.c)

    @property
    def pos(self):
        return self.c.pos

    def __repr__(self):
        return """HileUp("{position}", counts:{vals})""".format(
                position=self.c.pos, ref=self.c.reference_base,
                vals=dict(Counter(self.bases)))

def query_pos(AlignedSegment b):
    return qpos(b._delegate)

cdef class Config:
    cdef hile_config_t c

    def __init__(self, tags=(), track_read_names=True,
            track_reads = False,
            track_base_qualities=False, track_mapping_qualities=False,
            min_base_quality=10, min_mapping_quality=10, include_flags=0,
            exclude_flags=1796):
        self.c = hile_init_config()
        self.c.track_mapping_qualities = track_mapping_qualities
        self.c.track_base_qualities = track_base_qualities
        self.c.track_read_names = track_read_names
        self.c.track_reads = track_reads
        self.c.min_base_quality = min_base_quality
        self.c.min_mapping_quality = min_mapping_quality
        self.c.include_flags = include_flags
        self.c.exclude_flags = exclude_flags

        assert len(tags) <= 2, "must specify 2 or fewer tags"
        for i, t in enumerate(tags):
            assert len(t) == 2, "tags must have length 2"
            self.c.tags[2*i] = ord(t[0])
            self.c.tags[2*i+1] = ord(t[1])


def pileup(AlignmentFile bam, str chrom, int position, Config cfg, ignore_base='@'):
    cdef HileUp hile = HileUp()
    cdef char ib = ord(ignore_base[0])
    hile.c = hileup(bam.htsfile, bam.header.ptr, bam.index, chrom.encode(),
            position, &cfg.c, ib)
    return hile

def example():
    cdef AlignmentFile samfile = AlignmentFile("tests/three.bam", "rb")

    h = pileup(samfile, "1", 1585270, Config())
    print(h.bases)
    print(h.read_names)
    print(h.bqs)
    print(h.deletions)
    print(h.insertions)


if __name__ == "__main__":
    example()
