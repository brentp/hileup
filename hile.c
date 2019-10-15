#include "hile.h"
#include "khash.h"

KHASH_SET_INIT_STR(strset)

static inline void hile_realloc(hile *h, hile_config_t *cfg) {
    if(h->n < h->cap) { return; }
    h->cap = (h->cap == 0)? 2: (2 * h->cap);
    h->bases = realloc(h->bases, sizeof(hile_basestrand_t) * (h->cap));
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
#ifdef HILE_VERBOSE
        fprintf(stderr, "[hileup] tag %c%c not found\n", tag[0], tag[1]);
#endif
    }
    if(!added) {
	    if(append) {
		int oldlen = strlen(h->tags[h->n-1]);
		h->tags[h->n-1] = realloc(h->tags[h->n-1], sizeof(char) * (3 + oldlen));
		strcpy(h->tags[h->n-1] + oldlen, "/.");
	    } else {
		h->tags[h->n-1] = malloc(sizeof(char) * 2);
		strcpy(h->tags[h->n-1], ".");
	    }
    }
}

hile_config_t hile_init_config(void) {
    hile_config_t c;
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


void fill(hile *h, bam1_t *b, int position, hile_config_t *cfg) {
    if(b->core.qual < cfg->min_mapping_quality){ return; }
    if((cfg->include_flags != 0) && ((cfg->include_flags & b->core.flag) != cfg->include_flags)) { return; }
    if((cfg->exclude_flags & b->core.flag) != 0) { return; }

    int r_off = b->core.pos;
    int q_off = 0;
    bool skip_last = false;

    uint32_t *cig = bam_get_cigar(b);
    uint32_t i;
    for(i=0; i < b->core.n_cigar; i++){
        skip_last = false;
        uint32_t element = cig[i];
        uint32_t op = element & BAM_CIGAR_MASK;
        uint32_t oplen = element >> BAM_CIGAR_SHIFT;
        bool consumes_ref = (bam_cigar_type(op) & 2) != 0;
        bool consumes_query = (bam_cigar_type(op) & 1) != 0;
        // insertions and deletions get assinged to previous entry.
        if(r_off == position + 1 && !skip_last) {
            if(op == BAM_CDEL) {
              h->deletions = realloc(h->deletions, sizeof(hile_deletion_t) * (h->n_deletions+1));
              h->deletions[h->n_deletions].index = h->n - 1;
              h->deletions[h->n_deletions].length = oplen;
              h->n_deletions++;
            } else if (op == BAM_CINS) {
              h->insertions = realloc(h->insertions, sizeof(hile_insertion_t) * (h->n_insertions+1));
              h->insertions[h->n_insertions].index = h->n-1;
              h->insertions[h->n_insertions].length = oplen;
              h->insertions[h->n_insertions].sequence = malloc((oplen + 1) * sizeof(char));
              for(int k=0; k < oplen; k++){
                  h->insertions[h->n_insertions].sequence[k] = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), k + q_off)];
              }
              h->insertions[h->n_insertions].sequence[oplen] = '\0';
              h->n_insertions++;
            }
        }

        if(r_off > position){ break; }
        if(consumes_query) { q_off += oplen; }
        if(consumes_ref) { r_off += oplen; }
        if(r_off < position){ continue; }

        if(!(consumes_query && consumes_ref)){ continue; }

        int over = r_off - position;
        if(over > q_off) { break;}
        if(over < 0) { continue; }
        uint8_t bq = 0;
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
        hile_basestrand_t bs;
        bs.base = "=ACMGRSVTWYHKDBN"[bam_seqi(seq, i)];
        bs.reverse_strand = (b->core.flag & BAM_FREVERSE) ? 1: 0;
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

hile *hile_init(void) {
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
    int i;
    if(h->tags != NULL) {
        for(i=0; i < h->n; i++) {
            free(h->tags[i]);
        }
        free(h->tags);
    }
    if(h->read_names != NULL) {
        for(i=0; i < h->n; i++) {
            free(h->read_names[i]);
        }
        free(h->read_names);
    }
    if(h->mqs != NULL) { free(h->mqs); }
    if(h->bqs != NULL) { free(h->bqs); }
    if(h->deletions != NULL) { free(h->deletions); }
    if(h->insertions != NULL) {
        for(i=0;i<h->n_insertions;i++) {
           free(h->insertions[i].sequence);
	}
	free(h->insertions);
    }
    if(h->bases != NULL) { free(h->bases); }
    free(h);
}


hile *hileup(htsFile *htf, bam_hdr_t *hdr, hts_idx_t *idx, const char *chrom, int position, hile_config_t *cfg) {
  int tid = bam_name2id(hdr, chrom);
  if(tid == -1){
#ifdef HILE_VERBOSE
    fprintf(stderr, "[hile] unknown chromosome %s\n", chrom);
#endif
    return NULL;
  }

  // track overlapping reads.
  khash_t(strset) *seen;
  seen = kh_init(strset);
  khint_t k;
  int absent;

  hts_itr_t *itr = sam_itr_queryi(idx, tid, position, position + 1);
  if (itr == NULL) {
#ifdef HILE_VERBOSE
    fprintf(stderr, "[hileup] unable to access region %s:%d", chrom, position);
#endif
    return NULL;
  }
  int slen;
  bam1_t *b = bam_init1();
  hile *h = hile_init();
#define is_primary(b) (!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)))

  while((slen = sam_itr_next(htf, itr, b)) > 0){
     k = kh_get(strset, seen, bam_get_qname(b));
     if(k != kh_end(seen)){
        free((char *)kh_key(seen, k));
        kh_del(strset, seen, k);
        continue;
     }
     fill(h, b, position, cfg);
     if(!is_primary(b) || bam_endpos(b) <= b->core.mpos || b->core.mpos > position || b->core.tid != b->core.mtid || b->core.pos > b->core.mpos){
       continue;
     }
     char *qname = bam_get_qname(b);
     k = kh_put(strset, seen, qname, &absent);
     if (absent) { kh_key(seen, k) = strdup(qname); }
  }

  for(k=0; k < kh_end(seen); ++k){
      if(kh_exist(seen, k)) {
          free((char *)kh_key(seen, k));
      }
  }
  kh_destroy(strset, seen);

  bam_destroy1(b);
  hts_itr_destroy(itr);
  return h;
}

int example(void) {
    //htsFile *htf = hts_open("/data/human/hg002.cram", "rC");
    htsFile *htf = hts_open("tests/soft.bam", "rb");
    int start = 10080;
    bam_hdr_t *hdr = sam_hdr_read(htf);
    //hts_idx_t *idx = sam_index_load(htf, "/data/human/hg002.cram");
    hts_idx_t *idx = sam_index_load(htf, "tests/soft.bam");
    if(0 != hts_set_fai_filename(htf, "/data/human/g1k_v37_decoy.fa")) {
            fprintf(stderr, "cant set fai");
            return 2;
    }
    hile_config_t cfg = hile_init_config();
    cfg.track_base_qualities = true;
    cfg.track_mapping_qualities = true;
    cfg.track_read_names = true;
    cfg.min_base_quality = 10;
    cfg.min_mapping_quality = 10;

    hile* h = hileup(htf, hdr, idx, "1", start, &cfg);
    fprintf(stderr, "%s:%d ", "1", start);
    int i;
    for(i=0; i < h->n; i++){
        fprintf(stderr, "%c", (char)h->bases[i].base);
    }
    if(cfg.track_mapping_qualities) {
	    fprintf(stderr, " ");
	    for(i=0; i < h->n; i++){
		fprintf(stderr, "%c", (char)(h->bqs[i] + 33));
	    }
    }
    if(cfg.tags[0] != 0) {
	    fprintf(stderr, " ");
	    for(i=0; i < h->n; i++){
		fprintf(stderr, "%d:%s ", i, h->tags[i]);
	    }
    }
    fprintf(stderr, "\n");

    hile_destroy(h);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(htf);
    return 0;
}

int not_main(int argc, char *argv[]) {
    htsFile *htf = hts_open("tests/ins.bam", "rb");
    int start = 28588;
    bam_hdr_t *hdr = sam_hdr_read(htf);
    hts_idx_t *idx = sam_index_load(htf, "tests/ins.bam");
    hile_config_t cfg = hile_init_config();
    cfg.track_base_qualities = false;
    cfg.track_mapping_qualities = false;
    //cfg.track_read_names = true;
    cfg.tags[0] = 'C';
    cfg.tags[1] = 'B';

    for (int st = start - 4; st < start + 4; st++) {

	    hile* h = hileup(htf, hdr, idx, "1", st, &cfg);
	    fprintf(stderr, "%s:%d ", "1", st);
	    for(int i=0; i < h->n; i++){
		fprintf(stderr, "%c", (char)h->bases[i].base);
	    }
	    for(int i=0; i < h->n_insertions; i++){
		fprintf(stderr, " ins:%s ", h->insertions[i].sequence);
	    }
	    if(cfg.track_mapping_qualities) {
		    fprintf(stderr, " ");
		    for(int i=0; i < h->n; i++){
			fprintf(stderr, "%c", (char)(h->bqs[i] + 33));
		    }
	    }
	    /*
	    if(cfg.tags[0] != 0) {
		    fprintf(stderr, " ");
		    for(int i=0; i < h->n; i++){
			fprintf(stderr, "%d:%s ", i, h->tags[i]);
		    }
	    }
	    */
	    fprintf(stderr, "\n");

	    hile_destroy(h);
    }
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    hts_close(htf);
    return example();
}
