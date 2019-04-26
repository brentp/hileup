# distutils: sources = hile.c
# distutils: include_dirs = .
from pysam.libcalignmentfile cimport AlignmentFile
from libc.stdlib cimport malloc
cimport numpy as np
import numpy as np
cdef dt = np.dtype([('index', np.uint16, 1), ('len', np.uint32, 1)])

cdef class HileUp:
    cdef hile *c

    @property
    def bases(HileUp self):
        cdef int i
        cdef char* s = <char *>malloc((self.c.n+1) * sizeof(char))
        for i in range(0, self.c.n):
            #if show_strand and self.c.bases[i].reverse_strand == 1:
            #  s[i] = self.c.bases[i].base + 32
            #else:
            s[i] = self.c.bases[i].base
        s[self.c.n] = 0
        pys = s[:self.c.n]
        free(s)
        return pys


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
            result.append(self.c.read_names[i])
        return result

    @property
    def deletions(self):
        if self.c.n_deletions == 0: return []
        return np.asarray(<deletion_t[:self.c.n_deletions]>self.c.deletions, dtype=dt)

    @property
    def insertions(self):
        if self.c.n_insertions == 0: return []
        return np.asarray(<insertion_t[:self.c.n_insertions]>self.c.insertions, dtype=dt)

    @property
    def size(self):
        return self.c.n

    def __dealloc__(self):
        hile_destroy(self.c)

cdef class Config:
    cdef config_t c

    def __init__(self, tags=(), track_read_names=True,
            track_base_qualities=False, track_mapping_qualities=False,
            min_base_quality=10, min_mapping_quality=10, include_flags=0,
            exclude_flags=1796):
        self.c = hile_init_config()
        self.c.track_mapping_qualities = track_mapping_qualities
        self.c.track_base_qualities = track_base_qualities
        self.c.track_read_names = track_read_names
        self.c.min_base_quality = min_base_quality
        self.c.min_mapping_quality = min_mapping_quality
        self.c.include_flags = include_flags
        self.c.exclude_flags = exclude_flags

        assert len(tags) <= 2, "must specify 2 or fewer tags"
        for i, t in enumerate(tags):
            assert len(t) == 2, "tags must have length 2"
            self.c.tags[2*i] = ord(t[0])
            self.c.tags[2*i+1] = ord(t[1])


def pileup(AlignmentFile bam, str chrom, int position, Config cfg):
    cdef HileUp hile = HileUp()
    hile.c = hileup(bam.htsfile, bam.header.ptr, bam.index, chrom.encode(),
            position, &cfg.c)
    return hile

def example():
    cdef AlignmentFile samfile = AlignmentFile("tests/three.bam", "rb")

    h = pileup(samfile, "1", 1585270)
    print(h.bases)
    print(h.read_names)
    print(h.bqs)
    print(h.deletions)
    print(h.insertions)


if __name__ == "__main__":
    example()
