// TACO
// cgtf.c
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>

#include "cgtf.h"

#define TRUE 1
#define FALSE 0
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define GTF_DELIM "\t\r\n"
#define GTF_ATTR_DELIM ";"
#define GTF_KEY_VALUE_DELIM "= "
#define GTF_VALUE_QUOTE '"'
#define GTF_MAX_FIELD_SIZE 256

#define LOCUS_FMT_STRING "L%zu\t%s\t%zu\t%zu\t%" PRIu64 "\t%zu\n"

gtf_t *
gtf_new() {
    int errnum;
    gtf_t *gtfp = malloc(sizeof *gtfp);
    if (gtfp == NULL) {
        errnum = errno;
        fprintf(stderr, "Error in gtf_new(): %s\n", strerror(errnum));
    }
    return gtfp;
}


void
gtf_free(gtf_t *gtfp) {
    if (gtfp->line != NULL)
        free(gtfp->line);
    if (gtfp->h != NULL)
        kh_destroy(gtf_attr_hash, gtfp->h);
    free(gtfp);
}


int
gtf_parse(gtf_t *gtfp, char *s, int copy_str) {
    char *linep, *line;
    int errnum;

    if (copy_str) {
        // copy input string
        line = linep = strdup(s);
        if (line == NULL) {
            errnum = errno;
            fprintf(stderr, "Error in gtf_parse_line(): %s\n", strerror(errnum));
            free(line);
            return -1;
        }
        // store pointer so string can be freed later
        gtfp->line = line;
    } else {
        // do not copy string - original string will be mangled by parse
        line = linep = s;
        gtfp->line = NULL;
    }
    // split by tab into fields
    gtfp->seqname = strsep(&linep, GTF_DELIM);
    gtfp->source = strsep(&linep, GTF_DELIM);
    gtfp->feature = strsep(&linep, GTF_DELIM);
    gtfp->start = atoi(strsep(&linep, GTF_DELIM)) - 1;
    gtfp->end = atoi(strsep(&linep, GTF_DELIM));
    gtfp->score = strsep(&linep, GTF_DELIM);
    gtfp->strand = strsep(&linep, GTF_DELIM);
    gtfp->frame = strsep(&linep, GTF_DELIM);
    gtfp->attrs = strsep(&linep, GTF_DELIM);
    gtfp->h = NULL;
    return 0;
}


int
gtf_hash_attrs(gtf_t *gtfp) {
    char *attrp, *attr, *fieldp, *key, *value, *quote;
    khash_t(gtf_attr_hash) *h;
    unsigned k;
    int ret;

    // init attributes hash table
	h = kh_init(gtf_attr_hash);

    // split attributes by semicolon
    attrp = gtfp->attrs;
    while ((attr = strsep(&attrp, GTF_ATTR_DELIM)) != NULL) {
        if (strlen(attr) == 0)
            break;
        // scan past whitespace
        fieldp = attr;
        while (isspace(*fieldp))
            fieldp++;
        // split key and value by whitespace
        key = strsep(&fieldp, GTF_KEY_VALUE_DELIM);
        value = strsep(&fieldp, GTF_KEY_VALUE_DELIM);
        // strip quotes
        value++;
        if ((quote = strchr(value, GTF_VALUE_QUOTE)) == NULL) {
            fprintf(stderr, "Error in gtf_hash_attrs(): value missing quotes\n");
            return -1;
        }
        *quote = '\0';
        // insert attribute into hash table
        k = kh_put(gtf_attr_hash, h, key, &ret);
        kh_val(h, k) = value;
        if (!ret) kh_del(gtf_attr_hash, h, k);
    }
    gtfp->h = h;
    return 0;
}

const char *
gtf_get_attr(gtf_t *gtfp, char *a) {
    khint_t k;
    int is_missing;
    // lazy load attributes since this is expensive
    if (gtfp->h == NULL)
        gtf_hash_attrs(gtfp);
    // lookup key
    k = kh_get(gtf_attr_hash, gtfp->h, a);
    is_missing = (k == kh_end(gtfp->h));
    if (is_missing)
        return NULL;
    if (!kh_exist(gtfp->h, k))
        return NULL;
    return kh_value(gtfp->h, k);
}

int
gtf_fprintf(gtf_t *gtfp, FILE *fp) {
    khint_t k;
    khash_t(gtf_attr_hash) *h;
    const char *key, *value;
    char sep = '\t';

    fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t",
            gtfp->seqname,
            gtfp->source,
            gtfp->feature,
            gtfp->start,
            gtfp->end,
            gtfp->score,
            gtfp->strand,
            gtfp->frame);
    if (gtfp->h != NULL) {
        h = gtfp->h;
        for (k = kh_begin(h); k != kh_end(h); ++k) {
            if (kh_exist(h, k)) {
                key = kh_key(h, k);
                value = kh_value(h, k);
                fprintf(fp, "%c%s \"%s\";", sep, key, value);
                sep = ' ';
            }
        }
    } else {
        fprintf(fp, "%c%s", sep, gtfp->attrs);
    }
    fprintf(fp, "\n");
    return 0;
}


int
window_overlap(const char *seqname1, size_t start1, size_t end1,
               const char *seqname2, size_t start2, size_t end2) {
    if (strcmp(seqname1, seqname2) != 0)
        return FALSE;
    return ((start1 < end2) && (start2 < end1));
}


int
gtf_index_loci(char *gtf_file, char *output_file) {
    char *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;

    FILE *infp = fopen(gtf_file, "r");
    FILE *outfp = fopen(output_file, "w");
    gtf_t gtf;
    char seqname[GTF_MAX_FIELD_SIZE] = {0};
    size_t start = 0;
    size_t end = 0;
    size_t numlines = 0;
    size_t locus_id = 1;
    uint64_t filepos = 0;

    // get first line to initialize window
    linelen = getline(&line, &linecap, infp);
    if (linelen == 0) {
        fprintf(stderr, "build_locus_index(): empty input file");
        fclose(infp);
        fclose(outfp);
        return -1;
    }
    gtf_parse(&gtf, line, FALSE);

    // initialize window
    strncpy(seqname, gtf.seqname, GTF_MAX_FIELD_SIZE);
    seqname[sizeof(seqname) -1] = '\0';
    start = gtf.start;
    end = gtf.end;
    numlines = 1;

    while ((linelen = getline(&line, &linecap, infp)) > 0) {
        gtf_parse(&gtf, line, FALSE);
        // check if feature is outside current window
        if (!window_overlap(seqname, start, end, gtf.seqname,
                            gtf.start, gtf.end)) {
            // print current locus
            fprintf(outfp, LOCUS_FMT_STRING, locus_id, seqname, start, end,
                    filepos, numlines);
            // reset window
            strncpy(seqname, gtf.seqname, GTF_MAX_FIELD_SIZE);
            seqname[sizeof(seqname) -1] = '\0';
            start = gtf.start;
            end = gtf.end;
            filepos = ((uint64_t) ftello(infp)) - linelen;
            numlines = 1;
            locus_id += 1;
        } else {
            // expand window
            if (gtf.end > end)
                end = gtf.end;
            numlines += 1;
        }
    }

    // print final window
    if (numlines > 0) {
        // print current locus
        fprintf(outfp, LOCUS_FMT_STRING, locus_id, seqname, start, end,
                filepos, numlines);
    }

    if (line != NULL)
        free(line);
    fclose(infp);
    fclose(outfp);
    return 0;
}


int
main() {
    // build_locus_index("/Users/mkiyer/Documents/lab/projects/taco/ccle/testaggregate/transfrags.gtf", "loci.txt");
    //
    // gtf_t *gtfp;
    //
    // // allocate gtf struct
    // if ((gtfp = gtf_new()) == NULL)
    //     return 1;
    //
    // gtf_parse_line(gtfp, "chr1\ttest\ttranscript\t100\t200\t0\t+\t.\texpr \"2.4953546671\"; transcript_id \"42.T4\"; ref \"0\"; sample_id \"42\";\n");
    //
    // printf("printout\n");
    // gtf_fprintf(gtfp, stdout);
    //
    // printf("printout2\n");
    //
    // const char * a;
    // a = gtf_get_attr(gtfp, "expr");
    // printf("attr: %s\n", a);
    // a = gtf_get_attr(gtfp, "transcript_id");
    // printf("attr: %s\n", a);
    //
    // gtf_fprintf(gtfp, stdout);
    // gtf_free(gtfp);

    return 0;
}
