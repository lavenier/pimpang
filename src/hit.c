
#include "struct.h"

/*
//-- struct
typedef struct {
  int num_seq;
  int pos_query;
  int pos_db;
  int diag;
} hit_t;
*/

/*
//-- struct
typedef struct {
  int nb_hit;
  hit_t *hit;
} hits_t;
*/



//-- function
hits_t *create_hits()
{
  hits_t *H = (hits_t *) malloc(sizeof(hits_t));
  H->nb_hit    = 0;
  H->hit = (hit_t *) malloc (sizeof(hit_t) * MAX_HIT_DPU);
  return H;
}

//-- function
void free_hits(hits_t *H)
{
  free(H->hit);
  free(H);
}

void find_dpu_hits(hits_t *HITS, uint32_t *prot, index_t *I, uint32_t *tab_offset, uint64_t *oprot, uint8_t *iprot)
{
  uint32_t offset;
  uint32_t prot_len = tab_offset[0];
  int64_t  SR, SL, XL, XR;
  __m64 M0F = _m_from_int64 (0xF0F0F0F00F0F0F0F);
  __m64 M01 = _m_from_int64 (0x0101010101010101);
  __m64 LEFT, RIGHT, Q_LEFT, Q_RIGHT;


  HITS->nb_hit = 0;

  for (int iq=SIZE_NEIGHBOR; iq < prot_len-(SEED_SIZE+SIZE_NEIGHBOR); iq += 1)
    {
      offset  = tab_offset[iq];
      if (I->pos_seq[offset] == 0xFFFF) continue;
      Q_LEFT  = _mm_set_pi32(prot[iq-SIZE_NEIGHBOR],prot[iq-SIZE_NEIGHBOR]);
      Q_LEFT  = _mm_and_si64 (Q_LEFT,M0F);
      Q_RIGHT = _mm_set_pi32(prot[iq+SEED_SIZE],prot[iq+SEED_SIZE]);
      Q_RIGHT = _mm_and_si64 (Q_RIGHT,M0F);
      while (I->pos_seq[offset] != 0xFFFF)
	{
	  LEFT = _mm_set_pi32(I->left[offset],I->left[offset]);
	  LEFT = _mm_and_si64 (LEFT,M0F);
	  LEFT = _mm_cmpeq_pi8 (LEFT,Q_LEFT);
	  LEFT = _mm_and_si64 (LEFT,M01);
	  XL   = _m_to_int64(LEFT);
	  SL   = _mm_popcnt_u64(XL);

	  RIGHT = _mm_set_pi32(I->right[offset],I->right[offset]);
	  RIGHT = _mm_and_si64 (RIGHT,M0F);
	  RIGHT = _mm_cmpeq_pi8 (RIGHT,Q_RIGHT);
	  RIGHT = _mm_and_si64 (RIGHT,M01);
	  XR    = _m_to_int64(RIGHT);
	  SR    = _mm_popcnt_u64(XR);

	  if (SL+SR >= 7)
	    {
	      int num_seq = I->num_seq[offset];
	      int pos_db = I->pos_seq[offset];
	      int idb = (int) oprot[num_seq] + pos_db;
	      int idb8 = ((idb>>3)-4)<<3;
	      if (idb8 < 0)  idb8 = idb;
	      int delta = idb-idb8;
	      int cpt = 0;
	      int WD = 88;
	      uint8_t *buf = &iprot[idb8];
	      int k = iq - delta;
	      for (int i=0; i<WD; i++)
		{
		  if (k>=0)
		    if (buf[i]==(prot[k]&0xF)) cpt += 1;
		  k += 1;
		}
	      if (cpt > 25)
		{
		  if (HITS->nb_hit<MAX_HIT_DPU)
		    {
		      HITS->hit[HITS->nb_hit].num_seq =  num_seq;
		      HITS->hit[HITS->nb_hit].pos_query = iq;
		      HITS->hit[HITS->nb_hit].pos_db = pos_db;
		      HITS->hit[HITS->nb_hit].diag = HITS->hit[HITS->nb_hit].pos_query - HITS->hit[HITS->nb_hit].pos_db;
		      HITS->nb_hit += 1;
		    }
		}
	    }
	  offset += 1;
	}
    }
}


void host_find_rank_hits(param_t *P, syncbuf_t *SBQuery, syncbuf_t *SBHit, database_t *DB)
{
  query_t *Q;
  while ((Q = syncbuf_get(SBQuery)) != NULL)
    {
      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	Q->HITS[num_dpu] = create_hits();

      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	find_dpu_hits(Q->HITS[num_dpu], Q->idx_prot, DB->index[num_dpu], Q->offset, DB->oprot[num_dpu], DB->iprot[num_dpu]);
      
      syncbuf_insert(SBHit,Q);
    }
  if (syncbuf_end(SBQuery)) syncbuf_insert(SBHit,NULL);
}


void pim_find_rank_hits(struct dpu_set_t set, param_t *P, syncbuf_t *SBQuery, syncbuf_t *SBHit, database_t *DB)
{
  query_t *Q;
  while ((Q = syncbuf_get(SBQuery))!=NULL)
    {
      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	Q->HITS[num_dpu] = create_hits();

      // broadcast index and query to DPU
      xfer_query_to_dpu(set,Q);
      // run DPU
      DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
      // transfert results
      xfer_dpu_to_res(P, set, Q->HITS);

      // DEBUG
      /*
	if (P->pim_mode)
	{
	struct dpu_set_t dpu;
	DPU_FOREACH(set, dpu) { DPU_ASSERT(dpu_log_read(dpu, stdout)); } // DEBUG
	}
      */
      syncbuf_insert(SBHit,Q);
    }
  if (syncbuf_end(SBQuery)) syncbuf_insert(SBHit,NULL);
}

//-- function
void compute_hits(syncbuf_t *SBQuery, syncbuf_t *SBHit, param_t *P, database_t *DB, struct dpu_set_t *set)
{
  if (P->pim_mode)
    {
      int k = SBQuery->nb_cons;
      #pragma omp parallel for num_threads(k)
      for (int i=0; i<k; i++)
	pim_find_rank_hits(set[i],P,SBQuery,SBHit,DB);
    }
  else
    {
      int k = SBQuery->nb_cons;
      #pragma omp parallel for num_threads(k)
      for (int i=0; i<k; i++)
	host_find_rank_hits(P,SBQuery,SBHit,DB);

    }
}
