// ================================================================================================
// PANG
// ================================================================================================         
//      
// Find Protein Alignment with No Gap between 2 protein banks (query vs data_base)      
//      
// Run on PIM and Multi-Core architectures      
//      
// D. Lavenier - 2023
// ================================================================================================         
//      
// COMPILATION    
// -----------    
//      
// Makefile       
//      
// BLASTP comparison        
// -----------------        
// The equivalent blastp command to get full aligment format output is :      
// blastp -query <query.fasta> -subject <db.fasta> -ungapped -out <output> -comp_based_stats F -evalue <evalue> -seg yes        
// to get tabular output :  
// blastp -query <query.fasta> -subject <db.fasta> -ungapped -out <output> -comp_based_stats F -evalue <evalue> -seg yes -outfmt 6
//      


#include "struct.h"

#ifndef DPU_BINARY
#define DPU_BINARY "./dpu_pang"
#endif

#include <time.h>
#include <sys/time.h>
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

double duration(struct timeval begin)
{
  struct timeval end;
  gettimeofday(&end, 0);
  long seconds = end.tv_sec - begin.tv_sec;
  long microseconds = end.tv_usec - begin.tv_usec;
  double elapsed = seconds + microseconds*1e-6;
  return elapsed;
}

int main (int argc, char *argv[])
{
  omp_set_nested(1);

  param_t *Param = get_parameter(argc,argv);
  int8_t  *SubMat = init_matrix(Param);
  
  if (Param->verbose)
    {
      printf ("----\n");
      printf ("PANG: Protein Alignment with No Gap\n");
      printf ("----\n");
      printf ("Configuration\n");
      if (Param->pim_mode) printf ("  - PIM mode activated (%d ranks)\n",Param->nb_ranks);
      printf ("  - %d host threads\n",Param->nb_threads);
    }

  FILE *fdb;
  if ((fdb=fopen(Param->db_file_name,"r"))==NULL)
    {
      fprintf(stderr,"\nERROR: cannot open database file name (%s)\nexit...\n\n",Param->db_file_name);
      exit(0);
    }

  if (Param->verbose) printf ("Load database\n");
  struct timeval start;
  gettimeofday(&start, 0);
  db_info_t *db_info = get_db_info(Param->db_file_name); 
  database_t *DB = create_database(db_info->nb_seq,db_info->nb_aa);
  if (Param->verbose) printf ("  - %d proteins (%ld aa)\n",db_info->nb_seq, db_info->nb_aa);
  if (Param->verbose) printf ("  - %d proteins / dpu\n",DB->nb_seq_dpu);
  if (((db_info->nb_seq / NB_DPU_PER_RANK) + 1) > MAX_NB_SEQ_DPU)
    {
      fprintf(stderr,"\nERROR: too much sequence per DPU (max = %d)\nexit...\n\n",MAX_NB_SEQ_DPU);
      exit(0);
    }
  load_database(fdb,DB);
  fclose(fdb);
  if (Param->verbose) printf ("  - %.2f seconds\n",duration(start));

  
  FILE *fquery;
  if ((fquery=fopen(Param->query_file_name,"r"))==NULL)
    {
      fprintf(stderr,"\nERROR: cannot open query file name (%s)\nexit...\n\n",Param->query_file_name);
      exit(0);
    }

  if (Param->verbose) printf ("Index database\n");
  gettimeofday(&start, 0);
  index_database(Param,DB);
  if (Param->verbose) printf ("  - index size/dpu: %d (max=%d)\n",DB->size_db_dpu,MAX_NB_INDEX_LINE_DPU);
  if (Param->verbose) printf ("  - %.2f seconds\n",duration(start));

  syncbuf_t *SBQuery;
  syncbuf_t *SBHit;
  syncbuf_t *SBAlign;
  
  if (Param->pim_mode)
    {
      SBQuery = create_syncbuf("Query",64,Param->nb_ranks);
      SBHit   = create_syncbuf("Hit  ",256,Param->nb_threads);
    }
  else
    {
      int kq = Param->nb_threads / 3;
      kq = kq*2 + 1;
      int kh = Param->nb_threads - kq;
      if (kh == 0) kh = 1;
      SBQuery = create_syncbuf("Query",64,kq);
      SBHit   = create_syncbuf("Hit  ",256,kh);
    }
  SBAlign = create_syncbuf("Align",256,1);

  FILE *falign;
  if ((falign=fopen(Param->align_file_name,"w"))==NULL)
    {
      fprintf(stderr,"\nERROR: cannot open database file name (%s)\nexit...\n\n",Param->align_file_name);
      exit(0);
    }
  
  struct dpu_set_t set[Param->nb_ranks];
  if (Param->pim_mode)
    {
      // allocate DPUs
      for (int i=0; i<Param->nb_ranks; i++)
	DPU_ASSERT(dpu_alloc_ranks(1, NULL, &set[i]));
      // load DPU program
      for (int i=0; i<Param->nb_ranks; i++)
	DPU_ASSERT(dpu_load(set[i], DPU_BINARY, NULL));
      // transfer database to DPUs
      if (Param->verbose) printf ("Transfer database to DPUs\n");
      gettimeofday(&start,0);
      #pragma omp parallel for num_threads (Param->nb_ranks)
      for (int i=0; i<Param->nb_ranks; i++)
	xfer_db_to_dpu(Param, set[i], DB);
      if (Param->verbose) printf ("  - done\n");
      if (Param->verbose) printf ("  - %.2f seconds\n",duration(start));
    }

  gettimeofday(&start, 0);
  if (Param->verbose) printf ("Process queries\n");
  if (Param->format == 99) exit(0);
  
  #pragma omp parallel sections num_threads(4)
  {
    #pragma omp section
    get_queries(SBQuery,fquery,Param,DB);
    #pragma omp section
    compute_hits(SBQuery,SBHit,Param,DB,set);
    #pragma omp section
    compute_aligns(SBHit,SBAlign,Param,DB,SubMat);
    #pragma omp section
    save_aligns(SBAlign,Param,falign);
  }

  fclose(fquery);
  fclose(falign);
  if (Param->pim_mode)
    for (int i=0; i<Param->nb_ranks; i++) DPU_ASSERT(dpu_free(set[i]));

  if (Param->verbose) printf ("  - %.2f seconds\n",duration(start));
  
  if (Param->verbose)
    {
      printf ("Statistics\n");
      syncbuf_printstat(SBQuery);
      syncbuf_printstat(SBHit);
      syncbuf_printstat(SBAlign);
    }
  
  syncbuf_free(SBQuery);
  syncbuf_free(SBHit);
  syncbuf_free(SBAlign);

  //for (int i=0; i<Param->nb_ranks; i++) dpu_free(set[i]);
  //free_database(DB);
  
  return 0;
}
