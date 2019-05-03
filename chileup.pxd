from pysam.libchtslib cimport *

cdef extern from "stdbool.h":
    ctypedef bint bool

cdef extern from "hile.h" nogil:
    ctypedef struct hile_config_t:
        uint8_t min_mapping_quality
        uint8_t min_base_quality
        uint16_t exclude_flags
        uint16_t include_flags
        bool track_read_names
        bool track_base_qualities
        bool track_mapping_qualities
        char tags[4];

    ctypedef struct hile_deletion_t:
        uint16_t index
        uint32_t length
    ctypedef struct hile_insertion_t:
        uint16_t index
        # char *sequence
        uint32_t length
    ctypedef struct hile_basestrand_t:
        uint8_t reverse_strand
        uint8_t base

    ctypedef struct hile:
        uint32_t pos
        char reference_base
        hile_basestrand_t *bases
        uint16_t n
        uint16_t cap
        uint8_t *bqs
        uint8_t *mqs
        char **read_names
        char **tags;
        hile_insertion_t *insertions
        uint16_t n_insertions
        hile_deletion_t *deletions
        uint16_t n_deletions
    void hile_destroy(hile *h)
    hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, char *chrom, int position, hile_config_t *cfg)
    hile_config_t hile_init_config()
