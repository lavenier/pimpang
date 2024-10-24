
#include "struct.h"

/*
//-- struct
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
*/

/*
//-- struct
typedef struct {
  int nb_align;
  align_t **align;
} aligns_t;
*/

aligns_t *create_aligns(int nb_align)
{
  aligns_t *A = (aligns_t *) malloc (sizeof(aligns_t));
  A->nb_align = 0;
  A->align = (align_t **) malloc (sizeof(align_t*) * nb_align);
  return A;
}

//-- function
void free_aligns(aligns_t *A)
{
  for (int i=0; i<A->nb_align; i++) if (A->align[i] != NULL) free_align (A->align[i]);
  free(A->align);
  free(A);
}

//-- function
void free_align(align_t *A)
{
  free (A->match);
  free (A->query_subseq);
  free (A->db_subseq);
  free (A->db_com);
  free (A->query_com);
  free (A);
}

int cmp_align(const void *a, const void *b)
{
  align_t *pa = *(align_t**)a;
  align_t *pb = *(align_t**)b;
  return (pa->score < pb->score) - (pa->score > pb->score);
}

int diff_align(align_t *a1, align_t *a2)
{
  if (a1->score != a2->score) return 1;
  if (a1->query_start != a2->query_start) return 1;
  if (a1->db_start != a2->db_start) return 1;
  if (a1->len != a2->len) return 1;
  if (strcmp(a1->db_com,a2->db_com)!=0) return 1;
  return 0;
}

int cmphit(const void *a, const void *b)
{
  hit_t *A = (hit_t *)a;
  hit_t *B = (hit_t *)b;
  if (A->num_seq > B->num_seq) return -11;
  if (A->num_seq < B->num_seq) return 1;
  if (A->diag > B->diag) return -1;
  if (A->diag < B->diag) return 1;
  if (A->pos_query > B->pos_query) return -1;
  if (A->pos_query < B->pos_query) return 1;
  return 0;
}

int closehit(hit_t A, hit_t B)
{
  if (A.num_seq != B.num_seq) return 0;
  if (A.diag != B.diag) return 0;
  if (A.pos_query - B.pos_query > MAX_HIT_DIST) return 0;
  return 1;
}


void compute_full_align(align_t *A, protein_t *query_prot, protein_t *db_prot, int8_t *submat, param_t *P, uint64_t nb_aa)
{
  A->query_com = (char *) malloc(MAX_SIZE_DISPLAY_COM);
  A->db_com = (char *) malloc(MAX_SIZE_DISPLAY_COM);
  int i=0;
  for (i=0; i<MAX_SIZE_DISPLAY_COM-1; i++)
    {
      if (query_prot->com[i] == ' ') break;
      A->query_com[i] = query_prot->com[i];
    }
  while (i<MAX_SIZE_DISPLAY_COM-1) { A->query_com[i] = ' '; i += 1; }
  A->query_com[MAX_SIZE_DISPLAY_COM-1]='\0';
  for (i=0; i<MAX_SIZE_DISPLAY_COM-1; i++)
    {
      if (db_prot->com[i] == ' ') break;
      A->db_com[i] = db_prot->com[i];
    }
  while (i<MAX_SIZE_DISPLAY_COM-1) { A->db_com[i] = ' '; i += 1; }
  A->db_com[MAX_SIZE_DISPLAY_COM-1]='\0';
  // copy alignment subsequences from query and db
  A->query_subseq = (uint8_t *) malloc(sizeof(uint8_t)*A->len);
  A->db_subseq = (uint8_t *) malloc(sizeof(uint8_t)*A->len);
  for (int i=0; i<A->len; i++)
    {
      A->query_subseq[i] = query_prot->seq[A->query_start+i];
      A->db_subseq[i] = db_prot->seq[A->db_start+i];
    }
  // compute identity and positive
  // compute score
  A->nb_ident = 0;
  A->nb_pos = 0;
  A->score = 0;
  A->match = (uint8_t *) malloc(sizeof(uint8_t)*A->len);
  int8_t v;
  for (int i=0; i<A->len; i++)
    {
      v = submat[(query_prot->seq[A->query_start+i]<<5)+db_prot->seq[A->db_start+i]];
      A->score += v;
      if (query_prot->seq[A->query_start+i]==db_prot->seq[A->db_start+i])
	{
	  A->nb_ident += 1; A->nb_pos += 1; A->match[i] = '|';
	}
      else
	{
	  if (v >= 0)
	    {
	      A->nb_pos += 1; A->match[i] = ':';
	    }
	  else
	    {
	      A->match[i] = ' ';
	    }
	}
    }
  // compute evalue, bitscore
  A->evalue = score2evalue(P,query_prot->length,nb_aa,A->score);
  A->bitscore = bitscore(P,A->score);
}

align_t *compute_score_align(protein_t *Q, protein_t *D, uint16_t pos_query, uint16_t pos_db, int8_t * SubMat)
{
  int iq = (int) pos_query;
  int id = (int) pos_db;
  int score = 0;
  int max_score = -1000;
  int iql = pos_query;
  int idl = pos_db;
  int iqr = pos_query;

  align_t *A = (align_t *) malloc(sizeof(align_t));

  for (int i=0; i<SEED_SIZE; i++)
    {
      score += SubMat[(Q->seq[iq]<<5)+D->seq[id]];
      if (score >= max_score) { max_score = score; iqr = iq; }
      if (max_score - score > XDROP) break;
      iq++; id++;
    }
  // test is score seed > 0  [seed anchor may be wrong due to complexity encoding]
  if (max_score <= 0)
    {
      A->score = 0;
      A->len = 0;
      A->query_start = pos_query;
      A->db_start = pos_db;
      return A;
    }
  
  while ((iq<Q->length) && (id<D->length))
    {
      score += SubMat[(Q->seq[iq]<<5)+D->seq[id]];
      if (score >= max_score) { max_score = score; iqr = iq; }
      if (max_score - score > XDROP) break;
      iq++; id++;
    }
  iq = pos_query-1;
  id = pos_db-1;
  score = max_score;
  while ((iq>=0) && (id>=0))
    {
      score += SubMat[(Q->seq[iq]<<5)+D->seq[id]];
      if (score >= max_score) { max_score = score; iql = iq; idl = id; }
      if (max_score - score > XDROP) break;
      iq--; id--;
    }
  A->score = max_score;
  A->query_start = iql;
  A->db_start = idl;
  A->len = iqr-iql + 1;

  return A;
}

void dpu_align (param_t *P, hits_t *HH, aligns_t *T_AA, protein_t *QP, int min_score, db_prot_t *DB_PROT, int8_t *SubMat, uint64_t nb_aa)
{
  //int spy1 = HH->nb_hit;
  // suppress close hits
  if (HH->nb_hit != 0)
    {
      qsort (HH->hit,HH->nb_hit,sizeof(hit_t),cmphit);
      int p = 0;
      for (int i=1; i<HH->nb_hit; i++)
	{
	  if (closehit(HH->hit[p],HH->hit[i])==0)
	    {
	      HH->hit[p+1].num_seq = HH->hit[i].num_seq;
	      HH->hit[p+1].pos_query = HH->hit[i].pos_query;
	      HH->hit[p+1].pos_db = HH->hit[i].pos_db;
	      HH->hit[p+1].diag = HH->hit[i].diag;
	      p += 1;
	    }
	}
      HH->nb_hit = p+1;
    }
  //int spy2 = HH->nb_hit;
  
  // compute alignment with no gap
  T_AA->nb_align = 0;
  for (int j=0; j<HH->nb_hit; j++)
    {
      int k = HH->hit[j].num_seq; // k = protein number in DB_PROT
      protein_t *DP = DB_PROT->protein[k];
      // compute score from hit positions in QP & DP
      //if (local_score(QP, DP, HH->hit[j].pos_query, HH->hit[j].pos_db, SubMat))
	{
	  align_t *A = compute_score_align(QP, DP, HH->hit[j].pos_query, HH->hit[j].pos_db, SubMat);
	  if (A->score < min_score)
	    {
	      free(A);
	    }
	  else
	    {
	      compute_full_align(A, QP, DP, SubMat, P, nb_aa);
	      T_AA->align[T_AA->nb_align] = A;
	      T_AA->nb_align += 1;
	    }
	}
    }
  //printf ("%d %d %d\n",spy1,spy2,T_AA->nb_align);
}

void rank_align(syncbuf_t *SBHit, syncbuf_t *SBAlign, param_t *P, database_t *DB, int8_t *SubMat)
{
  query_t *Q;

  while ((Q = syncbuf_get(SBHit))!=NULL)
    {
      // create align data structure
      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	Q->ALIGNS[num_dpu] = create_aligns(Q->HITS[num_dpu]->nb_hit);
      
      // compute min score
      int min_score = evalue2score(P,Q->protein->length,DB->nb_aa) - 1;
      
      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	{
	  dpu_align(P, Q->HITS[num_dpu], Q->ALIGNS[num_dpu], Q->protein, min_score, DB->db_prot[num_dpu], SubMat, DB->nb_aa);
	}
      
      // calculate total number of alignments
      int nb_tot_align = 0;
      for (int i=0; i<NB_DPU_PER_RANK; i++)
	nb_tot_align += Q->ALIGNS[i]->nb_align;
      
      // allocate alignment memory array for the current query 
      Q->align = (align_t **) malloc (sizeof(align_t *)*nb_tot_align);
      Q->nb_align = 0;
      Q->max_nb_align = nb_tot_align;
      
      // gather all alignments
      for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu++)
	{
	  for (int i=0; i<Q->ALIGNS[num_dpu]->nb_align; i++)
	    {
	      Q->align[Q->nb_align] = Q->ALIGNS[num_dpu]->align[i];
	      Q->nb_align += 1;
	    }
	}
      
      if (Q->nb_align > 0)
	{
	  // sort alignments, best score first
	  qsort(Q->align,Q->nb_align,sizeof(align_t *),cmp_align);
	  
	  // eliminate identical alignments
	  int p = 0;
	  for (int i=1; i<Q->nb_align; i++)
	    {
	      if (diff_align(Q->align[p], Q->align[i]))
		{
		  Q->align[p+1] = Q->align[i];
		  p += 1;
		}
	    }
	  Q->nb_align = p+1;
	}
      syncbuf_insert(SBAlign,Q);
    }
  if (syncbuf_end(SBHit)) syncbuf_insert(SBAlign,NULL);
}


//-- function
void compute_aligns (syncbuf_t *SBHit, syncbuf_t *SBAlign, param_t *P, database_t *DB, int8_t *SubMat)
{
  int k = SBHit->nb_cons;
  #pragma omp parallel for num_threads(k)
  for (int i=0; i<k; i++)
    rank_align(SBHit,SBAlign,P,DB,SubMat);
}

