`hileup` is an early-stage version of a pileup engine.

It aims to provide an interface that is:
+ easy-to-use
+ fast from an interpreted language (like python).

It is currently targetted for accessing targetted sites (e.g. < 100K sites), rather
than sweeping across every site in the genome.

There is a version in nim, one in C, and a cython wrapper for the C in python.
The python version, which takes a pysam AlignmentFile object looks like:

```Python
import pysam
import chileup

bam = pysam.AlignmentFile("tests/three.bam", "rb")

config = chileup.Config(tags=[], track_read_names=True,
        track_base_qualities=True, track_mapping_qualities=True,
        exclude_flags=pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY | pysam.FDUP,
        min_base_quality=10, min_mapping_quality=10)


h = chileup.pileup(bam, "1", 1585270, config)

print(h.bases) # 'TT'
print(h.read_names) # [b'A00227:74:HCWC7DSXX:1:1269:13449:13855', b'A00227:74:HCWC7DSXX:1:2426:7157:15483']
print(h.bqs) # [37, 37]  numpy array that is a view into underlying data.
print(h.mqs) # [60, 60] numpy view.
print(h.deletions) # [(0, 8) (1, 8)] numpy view
print(h.insertions) # numpy view
print(h.tags) # copy.
```

To build `python setup.py build_ext -i`
To install `python setup.py install`

Because it minimizes operations in python, it is quite fast (for python).

**NOTE** that currently the strand information is unavailable from python.
