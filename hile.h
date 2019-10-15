#ifndef HILE_H
#define HILE_H

#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "stdint.h"
#include "stdio.h"
#include "string.h"
#include "stdbool.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * config_t determines which reads and bases are added to the
 * pileup. It should be created with `hile_init_config`
*/
typedef struct {
  uint8_t min_mapping_quality;
  uint8_t min_base_quality;
  uint16_t exclude_flags;
  uint16_t include_flags;
  bool track_read_names;
  bool track_base_qualities;
  bool track_mapping_qualities;
  // BAQ https://doi.org/10.1093/bioinformatics/btr076
  // uses sam_prob_realn in htslib
  //bool adjust_base_quality;
  char tags[4]; // track up to 2 tags.
} hile_config_t;

/*
 * hile_deletion_t holds the length of the event.
 * the `index` allows lookup into the hile struct
 * bases,read_names,etc that this event follows.
 */
typedef struct {
  uint32_t index;
  uint32_t length;
} hile_deletion_t;

/*
 * hile_insertion_t holds the length of the event.
 * the `index` allows lookup into the hile struct
 * bases,read_names,etc that this event follows.
 */
typedef struct {
  uint32_t index;
  char *sequence;
  uint32_t length;
} hile_insertion_t;

/*
 * hile_basestrand_t holds the base (A/C/G/T) and
 * the strand in a single uint8
 */
typedef struct {
  uint8_t reverse_strand:1, base: 7;
} hile_basestrand_t;

/*
 * hile tracks the "pileup at a single
 * genomic location
 */
typedef struct {
  uint32_t pos;
  char reference_base;
  hile_basestrand_t *bases;
  uint32_t n;
  uint32_t cap;
  uint8_t *bqs;
  uint8_t *mqs;
  char **read_names;
  char **tags;
  hile_insertion_t *insertions;
  uint32_t n_insertions;
  hile_deletion_t *deletions;
  uint32_t n_deletions;
} hile;


/*
 * create a hile for the given genomic position and config
 * the user is responsible for freeing the returned hile
 * using `hile_destroy`
 */
hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, const char *chrom, int position, hile_config_t *cfg);
/* free all memory from the given hile struct, including the insertion sequences */
void hile_destroy(hile *h);
/* initialize a config struct with sane defaults */
hile_config_t hile_init_config(void);

#ifdef __cplusplus
}
#endif

#endif
