`hileup` is an early-stage version of a pileup engine.

It aims to provide an interface that is:
+ easy-to-use
+ fast from an interpreted language (like python).

It is currently targetted for accessing targetted sites (e.g. < 100K sites), rather
than sweeping across every site in the genome.

There is a version in nim, one in C, and a cython wrapper for the C in python.
The python version, which takes a pysam AlignmentFile object looks like:

```Python

from pysam import AlignmentFile
from hileup import pileup
# this
samfile = AlignmentFile("tests/three.bam", "rb")

h = pileup(samfile, "1", 1585270)
print(h.bases) # 'TTT'
print(h.read_names) # qname from the alignment in order of bases
print(h.bqs) # base-qualities encoded as a string (subtract 33 from each char to get qual)
print(h.deletions) # a list of tuples indicating the index in the h.bases string and the length.
print(h.insertions)
```

The `nim` implementation is fairly optimized, the C implementation does not
correct for read-overlaps.

But, because it minimizes operations in python, it is quite fast (for python).
