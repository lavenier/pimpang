
#include <stdio.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <seqread.h>
#include <built_ins.h>
#include <barrier.h>

#include "constant.h"

typedef struct {
  int num_seq;
  int pos_query;
  int pos_db;
  int diag;
} hit_t;

// index storage [in MRAM]
__mram_noinit uint32_t   I_left[MAX_NB_INDEX_LINE_DPU];        // 16 MB
__mram_noinit uint32_t   I_right[MAX_NB_INDEX_LINE_DPU];       // 16 MB
__mram_noinit uint16_t   I_num_seq[MAX_NB_INDEX_LINE_DPU];     // 8 MB
__mram_noinit uint16_t   I_pos_seq[MAX_NB_INDEX_LINE_DPU];     // 8 MB   

// data base
__mram_noinit uint64_t  oprot[MAX_NB_SEQ_DPU];    // 40 KB
__mram_noinit uint8_t    iprot[MAX_SIZE_DB_DPU];   // 8 MB

// hit results
__mram_noinit hit_t hit[MAX_HIT_DPU];

// pre computed offset
__host uint32_t tab_offset[MAX_QUERY_LEN_DPU];      // 4 KB WRAM  (1K x 4)
// encoded query sequence
__host uint32_t prot[MAX_QUERY_LEN_DPU];            // 4 KB WRAM  (1K x 4)    

// hit
__host int hit_nb; 

#define SIZE_CACHE_INDEX 32

MUTEX_INIT(hit_mutex);
BARRIER_INIT(start_barrier, NR_TASKLETS);

int main()
{
  // cache realted to index I stored in MRAM
  __dma_aligned uint32_t IC_left[SIZE_CACHE_INDEX];
  __dma_aligned uint32_t IC_right[SIZE_CACHE_INDEX];
  __dma_aligned uint16_t IC_num_seq[SIZE_CACHE_INDEX];
  __dma_aligned uint16_t IC_pos_seq[SIZE_CACHE_INDEX];

  // hit buffer
  __dma_aligned int hitbuf[4];

  // offset buffer
  __dma_aligned uint64_t offbuf[1];

  // sequence database buffet
  __dma_aligned uint8_t buf[128];
  
  uint32_t LEFT1, Q_LEFT1, RIGHT1, Q_RIGHT1, LEFT2, Q_LEFT2, RIGHT2, Q_RIGHT2;
  uint32_t L1, L2, R1, R2, SR, SL;
  uint32_t prot_len = tab_offset[0];

  if (me()==0) hit_nb = 0;
  barrier_wait(&start_barrier);
    
  for (int iq=SIZE_NEIGHBOR+me(); iq < prot_len-(SEED_SIZE+SIZE_NEIGHBOR); iq += NR_TASKLETS)
    {
      uint32_t offset  = tab_offset[iq];
      Q_LEFT1  = prot[iq-SIZE_NEIGHBOR] & 0x0F0F0F0F;
      Q_LEFT2  = prot[iq-SIZE_NEIGHBOR] & 0xF0F0F0F0;
      Q_RIGHT1 = prot[iq+SEED_SIZE] & 0x0F0F0F0F;
      Q_RIGHT2 = prot[iq+SEED_SIZE] & 0xF0F0F0F0;
      mram_read(&I_left[offset],    IC_left,    SIZE_CACHE_INDEX*sizeof(uint32_t));
      mram_read(&I_right[offset],   IC_right,   SIZE_CACHE_INDEX*sizeof(uint32_t));
      mram_read(&I_num_seq[offset], IC_num_seq, SIZE_CACHE_INDEX*sizeof(uint16_t));
      mram_read(&I_pos_seq[offset], IC_pos_seq, SIZE_CACHE_INDEX*sizeof(uint16_t));
      uint32_t p = 0;
      while (IC_pos_seq[p] != 0xFFFF)
	{
	  LEFT1 = IC_left[p] & 0x0F0F0F0F;
	  __builtin_cmpb4_rrr(L1,Q_LEFT1,LEFT1);
	  // simulate popcount
	  L1 = L1 + (L1>>16);
	  L1 = L1 + (L1>>8);
	  L1 = L1 & 0xF;
	  LEFT2 = IC_left[p] & 0xF0F0F0F0;
	  __builtin_cmpb4_rrr(L2,Q_LEFT2,LEFT2);
	  // simulate popcount
	  L2 = L2 + (L2>>16);
	  L2 = L2 + (L2>>8);
	  L2 = L2 & 0xF;
	  SL = L1 + L2;

	  RIGHT1 = IC_right[p] & 0x0F0F0F0F;
	  __builtin_cmpb4_rrr(R1,Q_RIGHT1,RIGHT1);
	  // simulate popcount
	  R1 = R1 + (R1>>16);
	  R1 = R1 + (R1>>8);
	  R1 = R1 & 0xF;
	  RIGHT2 = IC_right[p] & 0xF0F0F0F0;
	  __builtin_cmpb4_rrr(R2,Q_RIGHT2,RIGHT2);
	  // simulate popcount
	  R2 = R2 + (R2>>16);
	  R2 = R2 + (R2>>8);
	  R2 = R2 & 0xF;
	  SR = R1 + R2;

	  if (SL + SR >= 7)
            {
	      int num_seq = IC_num_seq[p];
	      int pos_db = IC_pos_seq[p];
	      mram_read(&oprot[num_seq], offbuf, sizeof(uint64_t));
	      uint64_t cidb = offbuf[0];
	      int idb = (int) cidb;
	      idb += pos_db;
	      int idb8 = ((idb>>3)-4)<<3;
	      if (idb8 < 0) idb8 = idb;
	      int delta = idb-idb8;
	      int cpt = 0;
	      int WD = 88;
	      mram_read(&iprot[idb8], buf, 128);
	      int k = iq - delta;
	      for (int i=0; i<WD; i++)
		{
		  if (k>=0)
		    if (buf[i]==(prot[k]&0xF)) cpt += 1;
		  k += 1;
		}
	      if (cpt > 25)
		{
		  // start mutex section
		  mutex_lock(hit_mutex);
		  if (hit_nb < MAX_HIT_DPU)
		    {
		      // typedef struct { int num_seq; int pos_query; int pos_db; int diag; } hit_t;
		      hitbuf[0] = num_seq;
		      hitbuf[1] = iq;
		      hitbuf[2] = pos_db;
		      hitbuf[3] = iq - pos_db;
		      mram_write(hitbuf,&hit[hit_nb],16);
		      hit_nb += 1;
		    }
		  mutex_unlock(hit_mutex);
		  // end mutex section
		}
	    }
	  p += 1;
	  if (p == SIZE_CACHE_INDEX)
	    {
	      offset += SIZE_CACHE_INDEX;
	      mram_read(&I_left[offset],    IC_left,    SIZE_CACHE_INDEX*sizeof(uint32_t));
	      mram_read(&I_right[offset],   IC_right,   SIZE_CACHE_INDEX*sizeof(uint32_t));
	      mram_read(&I_num_seq[offset], IC_num_seq, SIZE_CACHE_INDEX*sizeof(uint16_t));
	      mram_read(&I_pos_seq[offset], IC_pos_seq, SIZE_CACHE_INDEX*sizeof(uint16_t));
	      p = 0;
	    }
	}
    }
}  



