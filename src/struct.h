
#define _GNU_SOURCE
#include <sched.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <immintrin.h>
#include <omp.h>

#include <semaphore.h>

#include <dpu.h>
#include <assert.h>
#include <dpu_log.h>

#include "constant.h"


// param.c
// -------

typedef struct {
  uint32_t  pim_mode;                  // 1: activated; 0: non activated
  double    evalue;                    // e-value
  char     *db_file_name;              // file data base name
  char     *query_file_name;           // file query name
  char     *align_file_name;           // file for storing alignments
  char     *submat_name;               // substitution matrix name
  float     lambda;                    // Lambda blast statistical parameter
  float     K;                         // K blast statistical parameter
  int       format;                    // output format
  int       verbose;                   // verbose mode
  int       nb_threads;                // number of host threads
  int       nb_ranks;                  // number of ranks used on the UPMEM memory
} param_t;


// matrix.c
// --------


// protein.c
// ---------

typedef struct {
  char *com;
  uint8_t *seq;
  uint32_t length;
} protein_t;

char AA[32] ;
int CODE_AA [256] ;

// hit.c
// -----

typedef struct {
  int num_seq;
  int pos_query;
  int pos_db;
  int diag;
} hit_t;

typedef struct {
  int nb_hit;
  hit_t *hit;
} hits_t;


// align.c
// -------

typedef struct {
  int        score;               // score of the alignment
  double     evalue;              // e-value of the alignment
  double     bitscore;            // bitscore of the alignment
  int        len;                 // length of the alignment
  int        query_start;         // starting position of the alignment on query sequence
  int        db_start;            // starting position of the alignment on database sequence
  int        nb_pos;              // number of positive matches
  int        nb_ident;            // number of identical matches
  char      *db_com;              // comment of the database sequence
  char      *query_com;           // comment of the query sequence
  uint8_t   *query_subseq;        // alignment subsequence from query sequence
  uint8_t   *db_subseq;           // alignment subsequence from database sequence
  uint8_t   *match;               // alignment
} align_t;

typedef struct {
  int nb_align;
  align_t **align;
} aligns_t;


// query.c
// -------

typedef struct {
  int num_query;                 // protein number
  uint32_t  *idx_prot;           // index code
  protein_t *protein;            // query protein sequence
  uint32_t  *offset;             // offset in the index
  uint8_t   *cplx;               // complexity zone
  uint32_t  nb_align;            // number of alignments
  uint32_t  max_nb_align;        // max_number of alignments
  align_t   **align;             // tabs of alignments
  hits_t    **HITS;              // tabs of hits - temporary array attached to the query
  aligns_t  **ALIGNS;            // tabs of alignments - temporary array attached to the query
} query_t;


// syncbuf.c
// ---------

typedef struct {
  sem_t mutex;
  sem_t empty;
  sem_t full;
  char *name;
  int nb_cons;
  int nb_end;
  int nb_elt;;
  int nb_empty;
  int nb_full;
  int size_buf;
  int index_write;
  int index_read;
  query_t **buf;
} syncbuf_t;


// index.c
// -------

int INDEX_CODE_AA [32] ;
typedef struct {
  uint32_t nb_db_prot;  // number of protein sequences indexed
  uint32_t *left;       // left neighborhood of the seed
  uint32_t *right;      // right neighborhood of the seed
  uint16_t *num_seq;    // sequence number 
  uint16_t *pos_seq;    // sequence position
} index_t;


// database.c
// ----------

typedef struct {
  uint32_t nb_aa;
  uint32_t nb_seq;
  protein_t **protein;
} db_prot_t;

typedef struct {
  uint32_t nb_seq;
  uint64_t nb_aa;
  uint32_t nb_seq_dpu;
  uint32_t size_db_dpu;
  uint32_t *seed_offset;
  db_prot_t **db_prot;
  index_t **index;
  uint64_t **oprot;
  uint8_t **iprot;
  uint32_t *size_iprot;
} database_t;

typedef struct {
  uint64_t nb_aa;
  uint32_t nb_seq;
} db_info_t;


// stat.c
// ------


// save.c
// ------


// dpu_xfer.c
// ----------


// param.c
// -------

param_t *get_parameter (int argc, char *argv[]);

// matrix.c
// --------

int8_t *init_matrix(param_t *P);
void free_matrix(int8_t *M);

// protein.c
// ---------

void free_protein(protein_t *P);
uint8_t *filter_complexity(protein_t *P);
protein_t *fasta_get_protein(FILE *fastafile);
void display_protein(protein_t *P);

// hit.c
// -----

hits_t *create_hits();
void free_hits(hits_t *H);
void compute_hits(syncbuf_t *SBQuery, syncbuf_t *SBHit, param_t *P, database_t *DB, struct dpu_set_t *set);

// align.c
// -------

void free_aligns(aligns_t *A);
void free_align(align_t *A);
void compute_aligns (syncbuf_t *SBHit, syncbuf_t *SBAlign, param_t *P, database_t *DB, int8_t *SubMat);

// query.c
// -------

void free_query(query_t *Q);
void get_queries(syncbuf_t *SB, FILE *fastafile, param_t *Param, database_t *DB);

// syncbuf.c
// ---------

syncbuf_t *create_syncbuf(char *name, int size_buffer, int nb_consumer);
void syncbuf_free (syncbuf_t *SB);
void syncbuf_insert(syncbuf_t *SB, query_t *Q);
void *syncbuf_get(syncbuf_t *SB);
int syncbuf_end(syncbuf_t *SB);
void syncbuf_printstat(syncbuf_t *SB);

// index.c
// -------

uint32_t *index_code_seq(uint8_t *seq, int len);
index_t *create_index(int size_db);
void free_index(index_t *I);
void index_database(param_t *P, database_t *DB);

// database.c
// ----------

db_prot_t *create_db_prot(uint32_t nb_seq);
void free_db_prot(db_prot_t *DP);
database_t *create_database(uint32_t nb_seq, uint64_t nb_aa);
void free_database(database_t *DB);
void load_database(FILE *fdb, database_t *DB);
db_info_t *get_db_info(char *file_db_name);

// stat.c
// ------

double bitscore(param_t *P, int32_t score);
double score2evalue(param_t *P, int query_len, uint64_t nb_aa_db, int score);
uint32_t evalue2score(param_t *P, int query_len, uint64_t nb_aa_db);

// save.c
// ------

void save_aligns(syncbuf_t *SBAlign, param_t *P, FILE *falign);

// dpu_xfer.c
// ----------

void xfer_db_to_dpu(param_t *P, struct dpu_set_t set, database_t *DB);
void xfer_query_to_dpu(struct dpu_set_t set, query_t *QRY);
void xfer_dpu_to_res(param_t *P, struct dpu_set_t set, hits_t **HITS);
