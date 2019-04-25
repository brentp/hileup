#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "stdint.h"
#include "stdio.h"
#include "string.h"
#include "stdbool.h"

typedef struct {
  uint8_t min_mapping_quality;
  uint8_t min_base_quality;
  uint16_t exclude_flags;
  uint16_t include_flags;
  bool track_read_names;
  bool track_base_qualities;
  bool track_mapping_qualities;
  char tags[4]; // track up to 2 tags.
} config_t;

typedef struct {
  uint16_t index;
  uint32_t length;
} deletion_t;

typedef struct {
  uint16_t index;
  //char *sequence;
  uint32_t length;
} insertion_t;

typedef struct {
  uint8_t reverse_strand:1, base: 7;
} basestrand_t;

typedef struct {
  uint32_t pos;
  char reference_base;
  basestrand_t *bases;
  uint16_t n;
  uint16_t cap;
  uint8_t *bqs;
  uint8_t *mqs;
  char **read_names;
  char **tags;
  insertion_t *insertions;
  uint16_t n_insertions;
  deletion_t *deletions;
  uint16_t n_deletions;
} hile;

void hile_destroy(hile *h);
hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, char *chrom, int position, config_t *cfg);
config_t hile_init_config();
