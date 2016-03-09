//
//  caggregate.c
//  TACO
//
//  Created by Balaji Pandian on 2/23/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
//  Need to Relay Safety Checks on sizes of strings, etc.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include <Python.h>

// Industry-standard FNV-32bit prime and offset
#define FNV_32_PRIME ((uint32_t)0x01000193)
#define FNV_32_OFFSET (2166136261U)

// Hash table size is 2^22 * sizeof(Hash_node) ~100mb
#define HASH_TABLE_SIZE (4194304)
#define MAX_GTF_LINE_SIZE (1024)
#define MAX_TRANSCRIPT_ID_SIZE (256)

typedef struct Hash_node {
    uint32_t num_exons;
    uint32_t length;
    char* id;
    char* transcript_id;
    struct Hash_node* next;
} Hash_node;

/* This is an industry-standard implementation of the FNV-32bit hash */
uint32_t fnv_32_str(char *str)
{
    unsigned char *s = (unsigned char *)str;
    uint32_t hval = FNV_32_OFFSET;

    while (*s) {
        hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
        hval ^= (uint32_t)*s++;
    }

    return hval % HASH_TABLE_SIZE;
}

bool str_in_hashtable(Hash_node* hash_table_root,char* str) {
    Hash_node bucket = hash_table_root[fnv_32_str(str)];

    if (bucket.transcript_id == NULL) {
        return false;
    }

    if (strcmp(bucket.transcript_id, str) == 0) {
        return true;
    }

    while (bucket.next != NULL) {
        if (strcmp(bucket.transcript_id, str) == 0) {
            return true;
        }

        bucket = *bucket.next;
    }

    return false;
}

void add_bucket_to_hashtable(Hash_node* hash_table_root,
                             char* transcript_string,
                             char* id_string) {
    char* allocated_transcript_string = strdup(transcript_string);
    char* allocated_id_string = strdup(id_string);
    uint32_t hash_value = fnv_32_str(transcript_string);

    if (hash_table_root[hash_value].transcript_id == '\0') {
        // add to initial bucket
        hash_table_root[hash_value].transcript_id = allocated_transcript_string;
        hash_table_root[hash_value].id = allocated_id_string;
        hash_table_root[hash_value].num_exons = 0;
        hash_table_root[hash_value].length = 0;
        hash_table_root[hash_value].next = NULL;
    } else {
        // add to linked list
        Hash_node* bucket = &hash_table_root[hash_value];
        while (bucket->next != NULL) {
            bucket = bucket->next;
        }

        Hash_node* new_hash_node = (Hash_node*)malloc(sizeof(Hash_node));
        new_hash_node->num_exons = 0;
        new_hash_node->length = 0;
        new_hash_node->next = NULL;
        new_hash_node->transcript_id = allocated_transcript_string;
        new_hash_node->id = allocated_id_string;

        bucket->next = new_hash_node;
    }
}

char* increment_node_and_return_id (Hash_node* hash_table_root, char* transcript_id, uint32_t length_increment) {
    Hash_node bucket = hash_table_root[fnv_32_str(transcript_id)];

    if (strcmp(bucket.transcript_id, transcript_id) == 0) {
        // Have the correct bucket.
    } else {
        while (strcmp(bucket.transcript_id, transcript_id) != 0) {
            bucket = *(bucket.next);
        }
    }

    bucket.num_exons += 1;
    bucket.length += length_increment;

    return bucket.id;
}

void free_all_hashtable_nodes(Hash_node* hash_table_root) {
    for (uint32_t i = 0; i < HASH_TABLE_SIZE; i++) {
        Hash_node* temp = &hash_table_root[i];

        if (temp->transcript_id != NULL) {
            free(temp->transcript_id);
            free(temp->id);
        }

        // Jump past first bucket -- allocated all at once in beginning
        temp = temp->next;

        while (temp != NULL) {
            Hash_node* temp_next = temp->next;
            free(temp->transcript_id);
            free(temp->id);
            free(temp);
            temp = temp_next;
        }
    }

    free(hash_table_root);
}

static PyObject* py_caggregate(PyObject* self, PyObject* args) {
    char* gtf_file;
    char* sample_id;
    char* gtf_expr_attr;
    char* output_file;
    char* stats_file;
    char* is_ref;

    if (!PyArg_ParseTuple(args, "ssssss", &gtf_file, &sample_id, &gtf_expr_attr, &output_file, &stats_file, &is_ref)) {
        return NULL;
    }

    size_t bufflen = MAX_GTF_LINE_SIZE;
    char* buffer = (char*)malloc(bufflen * sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    unsigned long long cur_t_id = 1;

    FILE* gtf_file_handler = fopen(gtf_file, "r");
    if (gtf_file_handler == NULL) {
        return NULL;
    }

    FILE* output_file_handler = fopen(output_file, "a");
    if (output_file_handler == NULL) {
        return NULL;
    }

    Hash_node* t_dict = (Hash_node*)calloc(HASH_TABLE_SIZE, sizeof(Hash_node));
    if (t_dict == NULL) {
        return NULL;
    }

    while (getline(&buffer, &bufflen, gtf_file_handler) != -1) {
        const char* seqname = strtok(buffer, "\t");
        const char* source = strtok(NULL, "\t");
        const char* feature = strtok(NULL, "\t");
        const char* start = strtok(NULL, "\t");
        const char* end = strtok(NULL, "\t");
        const char* score = strtok(NULL, "\t");
        const char* strand = strtok(NULL, "\t");
        const char* frame = strtok(NULL, "\t");
        const char* attribute = strtok(NULL, "\t");

        // Get Transcript ID from attributes
        char* transcript_id_start_index = strstr(attribute, "transcript_id \"");
        transcript_id_start_index += 15;
        char* end_of_transcript_id = strchr(transcript_id_start_index, '"');
        char transcript_id[MAX_TRANSCRIPT_ID_SIZE];
        strncpy(transcript_id, transcript_id_start_index, end_of_transcript_id - transcript_id_start_index);
        transcript_id[end_of_transcript_id - transcript_id_start_index] = '\0';

        if (strcmp(feature, "transcript") == 0) {
            char id_string[MAX_TRANSCRIPT_ID_SIZE];

            // Search for transcript_id in hashtable
            if (str_in_hashtable(t_dict, transcript_id)) {
                printf("GTF %s transcript_id %s is not unique", gtf_file, transcript_id);
                return NULL;
            } else {
                sprintf(id_string, "%s.T%llu", sample_id, cur_t_id);
                cur_t_id++;

                add_bucket_to_hashtable(t_dict, transcript_id, id_string);
            }

            if (strcmp(is_ref, "True") == 0) {
                fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texpr \"0.0\"; transcript_id \"%s\"; ref \"1\"; sample_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, id_string, sample_id);
            } else {

                // Get GTF_EXPR_ATTR (FPKM, for example) from attributes
                char* expr_start_index = strstr(attribute, gtf_expr_attr);

                // 2 for the space, then quotation mark
                expr_start_index += strlen(gtf_expr_attr) + 2;
                char* end_of_expr_index = strchr(expr_start_index, '"');
                char expr[64];
                strncpy(expr, expr_start_index, end_of_expr_index - expr_start_index);
                expr[end_of_expr_index - expr_start_index] = '\0';

                fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texpr \"%s\"; transcript_id \"%s\"; ref \"0\"; sample_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, expr, id_string, sample_id);
            }

        } else if (strcmp(feature, "exon") == 0) {
            char* id_string  = increment_node_and_return_id(t_dict, transcript_id, (uint32_t) (strtoul(end, NULL, 0) - strtoul(start, NULL, 0)));
            fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, id_string);
        }
    }

    fclose(gtf_file_handler);
    fclose(output_file_handler);
    free_all_hashtable_nodes(t_dict);
    free(buffer);

    // Return value -- returning "NULL" is an error
    return Py_BuildValue("i", 0);
}

// Mapping between python and c function names
static PyMethodDef CAggregateMethods[] = {
    {"caggregate", py_caggregate, METH_VARARGS},
    {NULL, NULL}
};

// Module initialisation routine
PyMODINIT_FUNC
initcaggregate(void)
{
    (void) Py_InitModule("caggregate", CAggregateMethods);
}
