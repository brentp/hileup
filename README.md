`hileup` is an early-stage version of a pileup engine.

It aims to provide an interface that is:
+ easy-to-use
+ fast from an interpreted language (like python).

It is currently targetted for accessing targetted sites (e.g. < 100K sites), rather
than sweeping across every site in the genome.

There is a version in nim, one in C, and a cython wrapper for the C in python.

## Python

The python version, which takes a pysam AlignmentFile object looks like:

```Python
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
```

To build `python setup.py build_ext -i`
To install `python setup.py install`

Because it minimizes operations in python, it is quite fast (for python).

**NOTE** that strand information is encoded by case for python (lower case == reverse strand).


## C

The C version should be transparent to anyone familier with [htslib](https://github.com/samtools/htslib)
The signature is:

```C
hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, char *chrom, int position, hile_config_t *cfg);
```

where `hile_config_t` is a simple struct that indicates min-mapping and base-qualities and whether to
track read-names, base-qualities, etc.


```C
    htsFile *htf = hts_open("tests/three.bam", "rb");
    int start = 1585270;
    bam_hdr_t *hdr = sam_hdr_read(htf);
    hts_idx_t *idx = sam_index_load(htf, "tests/three.bam");
    hile_config_t cfg = hile_init_config();
    cfg.track_base_qualities = true;
    cfg.track_mapping_qualities = true;
    cfg.track_read_names = true;
    // track the cell-barcode so we can get per-cell pileup!!
    cfg.tags[0] = 'C';
    cfg.tags[1] = 'B';

    hile* h = hileup(htf, hdr, idx, "1", start, &cfg);
    fprintf(stderr, "%s:%d ", "1", start);
    for(int i=0; i < h->n; i++){
        fprintf(stderr, "%c", (char)h->bases[i].base);
    }
    if(cfg.track_mapping_qualities) {
            fprintf(stderr, " ");
            for(int i=0; i < h->n; i++){
                fprintf(stderr, "%c", (char)(h->bqs[i] + 33));
            }
    }
    if(cfg.tags[0] != 0) {
            fprintf(stderr, " ");
            for(int i=0; i < h->n; i++){
                fprintf(stderr, "%d:%s ", i, h->tags[i]);
            }
    }
    fprintf(stderr, "\n");

    hile_destroy(h);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(htf);
```
