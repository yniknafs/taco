#include "klib/khash.h"

KHASH_MAP_INIT_STR(gtf_attr_hash, char*)

typedef struct gtf {
    char *line;
    char *seqname;
    char *source;
    char *feature;
    int start;
    int end;
    char *score;
    char *strand;
    char *frame;
    char *attrs;
    khash_t(gtf_attr_hash) *h;
} gtf_t;

int gtf_index_loci(char *gtf_file, char *output_file);
