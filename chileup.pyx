# distutils: sources = hile.c
# distutils: include_dirs = .
from pysam.libcalignmentfile cimport AlignmentFile
from libc.stdlib cimport malloc

cdef class HileUp:
    cdef hile *c

    @property
    def bases(self):
        cdef int i
        cdef char* s = <char *>malloc((self.c.n+1) * sizeof(char))
        for i in range(0, self.c.n):
            s[i] = self.c.bases[i].base
        s[self.c.n] = 0
        return s

    @property
    def bqs(self):
        cdef int i
        cdef char* s = <char *>malloc((self.c.n+1) * sizeof(char))
        for i in range(0, self.c.n):
            s[i] = (self.c.bqs[i] + 33)
        s[self.c.n] = 0
        return s

    @property
    def read_names(self):
        cdef int i
        result = []
        for i in range(0, self.c.n):
            result.append(self.c.read_names[i])
        return result

    @property
    def deletions(self):
        cdef int i
        cdef deletion_t ii
        result = []
        for i in range(0, self.c.n_deletions):
            ii = self.c.deletions[i]
            result.append((ii.index, ii.length))
        return result

    @property
    def insertions(self):
        cdef int i
        cdef insertion_t ii
        result = []
        for i in range(0, self.c.n_insertions):
            ii = self.c.insertions[i]
            result.append((ii.index, ii.length))
        return result

    @property
    def size(self):
        return self.c.n

    def __dealloc__(self):
        hile_destroy(self.c)

def pileup(AlignmentFile bam, str chrom, int position):
    cdef config_t cfg
    cfg.min_base_quality = 10
    cfg.min_mapping_quality = 10
    cfg.include_flags = 0
    cfg.exclude_flags = 0
    cfg.track_read_names = True
    cfg.track_base_qualities = True
    cfg.track_mapping_qualities = True
    cdef HileUp hile = HileUp()
    hile.c = hileup(bam.htsfile, bam.header.ptr, bam.index, chrom.encode(),
            position, &cfg)
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
