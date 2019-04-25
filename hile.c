#include "hile.h"

static inline void hile_realloc(hile *h, config_t *cfg) {
    if(h->n < h->cap) { return; }
    h->cap = (h->cap == 0)? 2: (2 * h->cap);
    h->bases = realloc(h->bases, sizeof(basestrand_t) * (h->cap));
    if(cfg->track_mapping_qualities) {
        h->mqs = realloc(h->mqs, sizeof(uint8_t) * (h->cap));
    }
    if(cfg->track_base_qualities) {
        h->bqs = realloc(h->bqs, sizeof(uint8_t) * (h->cap));
    }
    if(cfg->track_read_names) {
        h->read_names = realloc(h->read_names, sizeof(char *) * (h->cap));
    }
    if(cfg->tags[0] != 0) {
        h->tags = realloc(h->tags, sizeof(char *) * (h->cap));
    }
}

void hile_add_tag(hile *h, bam1_t *b, char tag[2], bool append) {
    uint8_t *t = bam_aux_get(b, tag);
    char *tval;
    bool added = false;
    if(t != NULL && (t[0] == 'H' || t[0] == 'Z' || t[0] == 'A')) {
        if(t[0] == 'H' || t[0] == 'Z') {
            tval = bam_aux2Z(t);
        } else {
            tval = bam_aux2Z(t);
        }
        if(tval != NULL) {
	    added = true;
	    if(append) {
		int oldlen = strlen(h->tags[h->n-1]);
		h->tags[h->n-1] = realloc(h->tags[h->n-1], sizeof(char) * (2 + strlen(tval) + oldlen));
		h->tags[h->n-1][oldlen] = '/';
		strcpy(h->tags[h->n-1] + oldlen + 1, tval);
	    } else {
		h->tags[h->n-1] = malloc(sizeof(char) * (1 + strlen(tval)));
		strcpy(h->tags[h->n-1], tval);
	    }
	}
    }
    if(t == NULL) {
        fprintf(stderr, "[hileup] tag %c%c not found\n", tag[0], tag[1]);
    }
    if(!added) {
	    if(append) {
		int oldlen = strlen(h->tags[h->n-1]);
		h->tags[h->n-1] = realloc(h->tags[h->n-1], sizeof(char) * (2 + oldlen));
		h->tags[h->n-1][oldlen] = '/';
		h->tags[h->n-1][oldlen+1] = '.';
	    } else {
		h->tags[h->n-1] = malloc(sizeof(char) * 2);
		strcpy(h->tags[h->n-1], ".");
	    }
    }
}

config_t hile_init_config() {
    config_t c;
    c.min_mapping_quality = 10;
    c.min_base_quality = 10;
    c.track_read_names = false;
    c.track_base_qualities = false;
    c.track_mapping_qualities = false;
    c.include_flags = 0;
    c.exclude_flags = BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY;
    c.tags[0] = 0;
    c.tags[1] = 0;
    c.tags[2] = 0;
    c.tags[3] = 0;
    return c;
}


void fill(hile *h, bam1_t *b, int position, config_t *cfg) {
    if(b->core.qual < cfg->min_mapping_quality){ return; }
    if((cfg->include_flags != 0) && ((cfg->include_flags & b->core.flag) != cfg->include_flags)) { return; }
    if((cfg->exclude_flags & b->core.flag) != 0) { return; }

    int r_off = b->core.pos;
    int q_off = 0;
    bool skip_last = false;

    uint32_t *cig = bam_get_cigar(b);
    for(uint32_t i=0; i < b->core.n_cigar; i++){
        skip_last = false;
	uint32_t element = cig[i];
	uint32_t op = element & BAM_CIGAR_MASK;
	uint32_t oplen = element >> BAM_CIGAR_SHIFT;
	bool consumes_ref = (bam_cigar_type(op) & 2) != 0;
	bool consumes_query = (bam_cigar_type(op) & 1) != 0;
        // insertions and deletions get assinged to previous entry.
        if(r_off == position + 1 && !skip_last) {
            if(op == BAM_CDEL) {
              h->deletions = realloc(h->deletions, sizeof(deletion_t) * (h->n_deletions+1));
              h->deletions[h->n_deletions].index = h->n - 1;
              h->deletions[h->n_deletions].length = oplen;
              h->n_deletions++;
            } else if (op == BAM_CINS) {
              h->insertions = realloc(h->insertions, sizeof(insertion_t) * (h->n_insertions+1));
              h->insertions[h->n_insertions].index = h->n-1;
              h->insertions[h->n_insertions].length = oplen;
              h->n_insertions++;
            }
        }

	if(r_off > position){ break; }
	if(consumes_query) { q_off += oplen; }
	if(consumes_ref) { r_off += oplen; }

	if(!(consumes_query && consumes_ref)){ continue; }

	int over = r_off - position;
	if(over > q_off) { break;}
	if(over < 0) { continue; }
	uint8_t bq = 0;
	// TODO handle overlapping mates.
	if(cfg->min_base_quality > 0 || cfg->track_base_qualities) {
	    bq = bam_get_qual(b)[q_off - over];
	    if (bq < cfg->min_base_quality) {
	       skip_last = true;
	       continue;
	    }
	}
	h->n++;
	hile_realloc(h, cfg);
        if(cfg->track_base_qualities) {
          h->bqs[h->n-1] = bq;
        }

        uint8_t *seq = bam_get_seq(b);
        int i = q_off - over;
        basestrand_t bs;
        bs.base = "=ACMGRSVTWYHKDBN"[(seq[i >> 1] >> ((~i & 1) >> 2) & 0xF)];
        bs.reverse_strand = b->core.flag & BAM_FREVERSE ? 1: 0;
        h->bases[h->n-1] = bs;

        if(cfg->track_read_names) {
            h->read_names[h->n-1] = malloc(sizeof(char) * (b->core.l_qname));
            strncpy(h->read_names[h->n-1], bam_get_qname(b), sizeof(char) * (b->core.l_qname));
        }
        if(cfg->track_mapping_qualities) {
            h->mqs[h->n - 1] = b->core.qual;
        }
	if(cfg->tags[0] != 0) {
	   char tag[2] = {cfg->tags[0], cfg->tags[1]};
           hile_add_tag(h, b, tag, false);
	   if(cfg->tags[2] != 0) {
	     tag[0] = cfg->tags[2], tag[1] = cfg->tags[3];
             hile_add_tag(h, b, tag, true);
	   }
	}
    }
}

hile *hile_init() {
  hile *h = malloc(sizeof(hile));
  h->read_names = NULL;
  h->bqs = NULL;
  h->tags = NULL;
  h->mqs = NULL;
  h->deletions = NULL;
  h->n_deletions = 0;
  h->insertions = NULL;
  h->n_insertions = 0;
  h->bases = NULL;
  h->pos = 0;
  h->n = 0;
  h->cap = 0;
  return h;
}


void hile_destroy(hile *h) {
    if(h == NULL) { return; }
    if(h->tags != NULL) {
        for(int i=0; i < h->n; i++) {
            free(h->tags[i]);
        }
        free(h->tags);
    }
    if(h->read_names != NULL) {
        for(int i=0; i < h->n; i++) {
            free(h->read_names[i]);
        }
        free(h->read_names);
    }
    if(h->mqs != NULL) { free(h->mqs); }
    if(h->bqs != NULL) { free(h->bqs); }
    if(h->deletions != NULL) { free(h->deletions); }
    if(h->insertions != NULL) { free(h->insertions); }
    if(h->bases != NULL) { free(h->bases); }
    free(h);
}

hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, char *chrom, int position, config_t *cfg) {
  char buffer[50];
  sprintf(buffer, "%s:%d-%d", chrom, position, position + 1);
  hts_itr_t *itr = sam_itr_querys(idx, hdr, buffer);
  if (itr == NULL) {
    fprintf(stderr, "[hileup] unable to access region %s", buffer);
    exit(1);
  }
  int slen;
  bam1_t *b = bam_init1();
  hile *h = hile_init();
  while((slen = sam_itr_next(htf, itr, b)) > 0){
     fill(h, b, position, cfg);
  }
  bam_destroy1(b);
  hts_itr_destroy(itr);
  return h;
}

int main() {
    htsFile *htf = hts_open("tests/three.bam", "rb");
    int start = 1585270;
    bam_hdr_t *hdr = sam_hdr_read(htf);
    hts_idx_t *idx = sam_index_load(htf, "tests/three.bam");
    config_t cfg = hile_init_config();
    cfg.track_base_qualities = true;
    cfg.track_mapping_qualities = true;
    cfg.track_read_names = true;
    cfg.tags[0] = 'M';
    cfg.tags[1] = 'D';
    cfg.tags[2] = 'F';
    cfg.tags[3] = 'G';

    hile* h = hileup(htf, hdr, idx, "1", start, &cfg);
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
		fprintf(stderr, "%s ", h->tags[i]);
	    }
    }
    fprintf(stderr, "\n");

    hile_destroy(h);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(htf);
    return 0;
}
