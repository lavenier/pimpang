
#include "struct.h"

//-- extern
int INDEX_CODE_AA [32] = {12, 11,  8,  7,  6, 10, 10, 13, 15,  1,  1, 10,  3,  2, 14,  9,  9,  5,  4,  1,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//-- function
uint32_t *index_code_seq(uint8_t *seq, int len)
{
  uint32_t *iseq = malloc(sizeof(uint32_t)*len);
  for (int i=0; i<len-SIZE_NEIGHBOR; i++)
    {
      iseq[i] = 0;
      int k = 0;
      for (int j=0; j<8; j++)
	{
	  iseq[i] = iseq[i] | (INDEX_CODE_AA[seq[i+j]]<<k);
	  k += 4;
	}
    }
  return iseq;
}

/*
//-- struct
typedef struct {
  uint32_t nb_db_prot;  // number of protein sequences indexed
  uint32_t *left;       // left neighborhood of the seed
  uint32_t *right;      // right neighborhood of the seed
  uint16_t *num_seq;    // sequence number 
  uint16_t *pos_seq;    // sequence position
} index_t;
*/

//-- function
index_t *create_index(int size_db)
{
  index_t *I  = (index_t *)  malloc (sizeof(index_t));
  I->left     = (uint32_t *) malloc (sizeof(uint32_t) * size_db);
  I->right    = (uint32_t *) malloc (sizeof(uint32_t) * size_db);
  I->num_seq  = (uint16_t *) malloc (sizeof(uint16_t) * size_db);
  I->pos_seq  = (uint16_t *) malloc (sizeof(uint16_t) * size_db);
  return I;
}


//-- function
void free_index(index_t *I)
{
  free(I->left);
  free(I->right);
  free(I->num_seq);
  free(I->pos_seq);
  free(I);
}

/*
             uint64_t *oprot                  uint8_t *iprot   (list of coded protein sequence - aligned on 8 multiple)
             -------------                    ----------------------------------------------------------------
            | offset P0   |                  |                         |                |              |
             -------------                    ----------------------------------------------------------------
            | offset P1   |                  ^                         ^                ^
             -------------                   |                         |                |
            | offset P2   |              offset P0                  offset P1        offset P2 
             -------------
            |             |

 */


void code_db_prot(database_t *DB, int num_dpu)
{
  DB->oprot[num_dpu] = (uint64_t *) malloc(sizeof(uint64_t)*MAX_NB_SEQ_DPU);
  uint32_t offset = 0;
  uint32_t len = 0;
  uint8_t *seq;
  for (int num_seq=0; num_seq < DB->db_prot[num_dpu]->nb_seq; num_seq++)
    {
      offset += len;
      DB->oprot[num_dpu][num_seq] = (uint64_t) offset;
      len = DB->db_prot[num_dpu]->protein[num_seq]->length;
      len = ((len>>3) + 1) << 3;
    }
  DB->size_iprot[num_dpu] = offset + len; 
  DB->iprot[num_dpu] = (uint8_t *) calloc(MAX_SIZE_DB_DPU,sizeof(uint8_t));
  for (int num_seq=0; num_seq < DB->db_prot[num_dpu]->nb_seq; num_seq++)
    {
      len = DB->db_prot[num_dpu]->protein[num_seq]->length;
      seq = DB->db_prot[num_dpu]->protein[num_seq]->seq;
      offset = (uint32_t) DB->oprot[num_dpu][num_seq];
      for (int i=0; i<len; i++)
	{
	  DB->iprot[num_dpu][offset+i] = INDEX_CODE_AA[seq[i]];
	}
    }
}


/*
 index structure
 ===============

           uint32_t                       uint16_t          uint16_t         uint32_t        uint32_t                     
           seed_offset                    pos_seq           num_seq          left            right
           -----------                   ----------        ---------       ----------       ----------
   0x0000 |     0     |     0x0000      |  0xFFFF  |      |         |     |          |     |          |
   0x0001 |     8     |     0x0001      |          |      |         |     |          |     |          |
          |           |     0x0002      |          |      |         |     |          |     |          |
          |           |     0x0003      |          |      |         |     |          |     |          |
          |           |     0x0004      |          |      |         |     |          |     |          |
          :           :     0x0005      |          |      |         |     |          |     |          |
          :           :     0x0006      |          |      |         |     |          |     |          |
          |           |     0x0007      |          |      |         |     |          |     |          |
   0xFFFF |           |     0x0008      |      35  |      |   28    |     |0x57A42B8F|     |0x6D634AB9|
           -----------      0x0009      |     183  |      |   72    |     |0x*BAC5622|     |0xBF56A321|
                            0x000A      |  0xFFFF  |      |         |     |          |     |          |
                                        :          :      :         :     :          :     :          :

*/


void index_db_prot(index_t *I, db_prot_t *D, uint32_t *seed_offset)
{
  I->nb_db_prot = D->nb_seq;
  uint32_t *off = (uint32_t *) malloc (sizeof(uint32_t)*SIZE_OFFSET_INDEX);
  for (int i=0; i<SIZE_OFFSET_INDEX; i++) off[i] = seed_offset[i];
  for (int i=0; i<SIZE_OFFSET_INDEX; i++) I->pos_seq[off[i]] = 0xFFFF;
  int seed_code;
  // fill position table
  for (int num_seq=0; num_seq < D->nb_seq; num_seq++)
    {
      uint8_t *seq = D->protein[num_seq]->seq;
      uint32_t len = D->protein[num_seq]->length;
      if (len > SIZE_NEIGHBOR+SEED_SIZE)
	{
	  uint32_t *iseq = index_code_seq(seq,len);
	  for (int pos_seq=SIZE_NEIGHBOR; pos_seq < len - (SIZE_NEIGHBOR+SEED_SIZE); pos_seq++)
	    {
	      seed_code = iseq[pos_seq]&SEED_MASK;
	      if (seed_code != 0)
		{
		  I->pos_seq[off[seed_code]] = pos_seq;
		  I->num_seq[off[seed_code]] = num_seq;
		  I->left[off[seed_code]] = iseq[pos_seq-SIZE_NEIGHBOR];
		  I->right[off[seed_code]] = iseq[pos_seq+SEED_SIZE];
		  off[seed_code] += 1;
		  I->pos_seq[off[seed_code]] = 0xFFFF;
		}
	    }
	  free(iseq);
	}
    }
  free(off);
}

//-- function
void index_database(param_t *P, database_t *DB)
{
  uint32_t *counter = (uint32_t *) malloc(sizeof(uint32_t)*SIZE_OFFSET_INDEX); // count number of seeds for a dpu
  uint32_t *max_cnt = (uint32_t *) malloc(sizeof(uint32_t)*SIZE_OFFSET_INDEX); // count max number of seeds for all dpu
  for (int i=0; i<SIZE_OFFSET_INDEX; i++) max_cnt[i] = 0;
  int seed_code;
  for (int num_dpu=0; num_dpu<NB_DPU_PER_RANK; num_dpu +=1)
    {
      for (int i=0; i<SIZE_OFFSET_INDEX; i++) counter[i] = 0; // reset counter for DPU num_dpu
      db_prot_t *D = DB->db_prot[num_dpu];                    // D = pointer to the set of proteins of DPU num_dpu
      for (int num_seq=0; num_seq<D->nb_seq; num_seq += 1)    // iterate on each sequence
	{
	  uint8_t *seq = D->protein[num_seq]->seq;
	  uint32_t len = D->protein[num_seq]->length;
	  if (len > SIZE_NEIGHBOR+SEED_SIZE)
	    {
	      uint32_t *iseq = index_code_seq(seq,len);           // encoded sequence
	      for (int pos_seq=SIZE_NEIGHBOR; pos_seq < len - (SIZE_NEIGHBOR+SEED_SIZE); pos_seq++)
		{
		  seed_code = iseq[pos_seq]&SEED_MASK;
		  if (seed_code != 0) counter[seed_code] += 1;    // count number of seeds (except for seed = 0)
		}
	      free (iseq);
	    }
	}
      for (int i=0; i<SIZE_OFFSET_INDEX; i++)                 // update max counter
	if (counter[i] > max_cnt[i]) max_cnt[i] = counter[i];
    }

  
  for (int i=0; i<SIZE_OFFSET_INDEX; i++) max_cnt[i] = ((max_cnt[i]>>3)+1)<<3;   // count must be a multiple of 8 (MRAM ==> WRAM DPU constraints)
  DB->size_db_dpu = 0;
  for (int i=0; i<SIZE_OFFSET_INDEX; i++) DB->size_db_dpu += max_cnt[i];  // size_db_dpu = number of lines of the index
  
  DB->seed_offset[0] = 0;
  for (int i=1; i<SIZE_OFFSET_INDEX; i++)
    {
      DB->seed_offset[i] = DB->seed_offset[i-1]+max_cnt[i-1];  // compute offset of each seed
    }
  free (counter);
  free (max_cnt);

  // create index for each DPU
  // the offset index is the same for all DPUs
  for (int i=0; i<NB_DPU_PER_RANK; i++) DB->index[i] = create_index(DB->size_db_dpu);
  #pragma omp parallel for num_threads(P->nb_threads)
  for (int i=0; i<NB_DPU_PER_RANK; i++)
    {
      index_db_prot(DB->index[i],DB->db_prot[i],DB->seed_offset);
    }
  #pragma omp parallel for num_threads(P->nb_threads)
  for (int i=0; i<NB_DPU_PER_RANK; i++)
    {
      code_db_prot(DB,i);
    }
}

