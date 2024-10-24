
#include "struct.h"

/*
//-- struct
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
*/


uint32_t *query_offset(uint32_t *idx_prot, uint32_t len, uint32_t *offset, uint8_t *cplx)
{
  uint32_t *OFFSET = (uint32_t *) malloc(sizeof(uint32_t)*len);
  OFFSET[0] = (uint32_t) len;
  for (int i=SIZE_NEIGHBOR; i<len-(SEED_SIZE+SIZE_NEIGHBOR); i++)
    {
      if (cplx[i] == 1)
	{
	  OFFSET[i] = 0;
	}
      else
	{
	  OFFSET[i] = offset[idx_prot[i]&SEED_MASK];
	}
    }
  return OFFSET;
}

query_t *create_query(protein_t *P, int num, database_t *DB)
{
  query_t *Q = (query_t *) malloc (sizeof(query_t));
  Q->num_query = num;
  Q->protein = P;
  Q->cplx = filter_complexity(P);
  Q->idx_prot = index_code_seq(P->seq, P->length);
  Q->offset = query_offset(Q->idx_prot, P->length, DB->seed_offset, Q->cplx);
  Q->HITS = (hits_t **) malloc(sizeof(hits_t *)*NB_DPU_PER_RANK);
  Q->ALIGNS = (aligns_t **) malloc(sizeof(aligns_t *)*NB_DPU_PER_RANK);
  return Q;
}

//-- function
void free_query(query_t *Q)
{
  free (Q->idx_prot);
  free_protein(Q->protein);
  free (Q->offset);
  free (Q->cplx);
  for (int i=0; i<NB_DPU_PER_RANK; i++) free_aligns(Q->ALIGNS[i]);
  free (Q->ALIGNS);
  for (int i=0; i<NB_DPU_PER_RANK; i++) free_hits(Q->HITS[i]);
  free (Q->HITS);
  free (Q->align);
  free (Q);
}

//-- function
void get_queries(syncbuf_t *SB, FILE *fastafile, param_t *Param, database_t *DB)
{
  protein_t *P;
  query_t *Q = NULL;
  int num_query = 0;
  while ((P = fasta_get_protein(fastafile)) != NULL)
    {
      if (P->length < MAX_QUERY_LEN_DPU)
	{
	  Q = create_query(P,num_query,DB);
	  syncbuf_insert(SB,Q);
	  num_query += 1;
	}
    }
  syncbuf_insert(SB,NULL);
 }
