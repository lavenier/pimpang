
#include "struct.h"

//-- function
void xfer_db_to_dpu(param_t *P, struct dpu_set_t set, database_t *DB)
{
  struct dpu_set_t dpu;
  uint32_t each_dpu;
  //printf ("DB->size_db_dpu : %d\n",DB->size_db_dpu);
  // transfer left neighborhood
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->index[each_dpu]->left));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"I_left",0,sizeof(uint32_t)*DB->size_db_dpu,DPU_XFER_DEFAULT));
  // transfer right neighborhood
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->index[each_dpu]->right));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"I_right",0,sizeof(uint32_t)*DB->size_db_dpu,DPU_XFER_DEFAULT));
  // transfer sequence number
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->index[each_dpu]->num_seq));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"I_num_seq",0,sizeof(uint16_t)*DB->size_db_dpu,DPU_XFER_DEFAULT));
  // transfer position in sequence
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->index[each_dpu]->pos_seq));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"I_pos_seq",0,sizeof(uint16_t)*DB->size_db_dpu,DPU_XFER_DEFAULT));
  // transfer iprot
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->iprot[each_dpu]));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"iprot",0,sizeof(uint8_t)*MAX_SIZE_DB_DPU,DPU_XFER_DEFAULT));
  // transfer oprot
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,DB->oprot[each_dpu]));
  DPU_ASSERT (dpu_push_xfer(set,DPU_XFER_TO_DPU,"oprot",0,sizeof(uint64_t)*MAX_NB_SEQ_DPU,DPU_XFER_DEFAULT));
}


//-- function
void xfer_query_to_dpu(struct dpu_set_t set, query_t *QRY)
{
  int xfer_size = ((QRY->protein->length>>3)+1)<<3;
  // broadcast encoded protein
  DPU_ASSERT(dpu_broadcast_to(set,"prot",0,QRY->idx_prot,sizeof(uint32_t)*xfer_size,DPU_XFER_DEFAULT));
  // broadcast offset
  DPU_ASSERT(dpu_broadcast_to(set,"tab_offset",0,QRY->offset,sizeof(uint32_t)*xfer_size,DPU_XFER_DEFAULT));
}


//-- function
void xfer_dpu_to_res(param_t *P, struct dpu_set_t set, hits_t **HITS)
{
  struct dpu_set_t dpu;
  uint32_t each_dpu;
  // number of hits computed per DPU
  uint32_t *dpu_hit_nb = (uint32_t *) malloc(sizeof(uint32_t)*NB_DPU_PER_RANK);
  // get number of hits per DPU
  DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,&dpu_hit_nb[each_dpu]));
  DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "hit_nb", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));

  // compute the maximum 
  int maxi=0;
  for (int i=0; i<NB_DPU_PER_RANK; i++)
    {
      HITS[i]->nb_hit = dpu_hit_nb[i];
      maxi = (dpu_hit_nb[i] > maxi) ? dpu_hit_nb[i] : maxi;
    }
  if (maxi > 0)
    {
      int size_xfer = ((maxi>>2)+1)<<2;
      if (size_xfer > MAX_HIT_DPU) size_xfer = MAX_HIT_DPU;
      // get hits
      DPU_FOREACH (set,dpu,each_dpu) DPU_ASSERT(dpu_prepare_xfer(dpu,HITS[each_dpu]->hit));
      DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "hit", 0, sizeof(hit_t)*size_xfer, DPU_XFER_DEFAULT));
    }
  free(dpu_hit_nb);
}


