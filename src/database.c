
#include "struct.h"

/*
//-- struct
typedef struct {
  uint32_t nb_aa;
  uint32_t nb_seq;
  protein_t **protein;
} db_prot_t;
*/

//-- function
db_prot_t *create_db_prot(uint32_t nb_seq)
{
  db_prot_t *DP = (db_prot_t *) malloc(sizeof(db_prot_t));
  DP->nb_aa = 0;
  DP->nb_seq = 0;
  DP->protein = (protein_t **) malloc(sizeof(protein_t *)*nb_seq);
  return DP;
}

//-- function
void free_db_prot(db_prot_t *DP)
{
  for (int i=0; i<DP->nb_seq; i++) free_protein(DP->protein[i]);
  free (DP);
}

/*
//-- struct
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
*/

//-- function
database_t *create_database(uint32_t nb_seq, uint64_t nb_aa)
{
  database_t *DB = (database_t *) malloc (sizeof(database_t));
  DB->nb_aa = nb_aa;
  DB->nb_seq = nb_seq;
  DB->nb_seq_dpu = (DB->nb_seq / NB_DPU_PER_RANK) + 1;
  DB->db_prot = (db_prot_t **) malloc (sizeof(db_prot_t *)*NB_DPU_PER_RANK);
  DB->index = (index_t **) malloc (sizeof(index_t *)*NB_DPU_PER_RANK);
  DB->seed_offset = (uint32_t *) malloc(sizeof(uint32_t)*SIZE_OFFSET_INDEX);
  DB->oprot = (uint64_t **) malloc(sizeof(uint64_t *)*NB_DPU_PER_RANK);
  DB->iprot = (uint8_t **) malloc (sizeof(uint8_t *)*NB_DPU_PER_RANK);
  DB->size_iprot = (uint32_t *) malloc (sizeof(uint32_t)*NB_DPU_PER_RANK);
  return DB;
}

//-- function
void free_database(database_t *DB)
{
  for (int i=0; i<NB_DPU_PER_RANK; i++) free_db_prot(DB->db_prot[i]);
  for (int i=0; i<NB_DPU_PER_RANK; i++) free_index(DB->index[i]);
  for (int i=0; i<NB_DPU_PER_RANK; i++) free(DB->iprot[i]);
  for (int i=0; i<NB_DPU_PER_RANK; i++) free(DB->oprot[i]);
  free (DB->db_prot);
  free (DB->index);
  free (DB->seed_offset);
  free (DB->oprot);
  free (DB->iprot);
  free (DB->size_iprot);
  free (DB);
}

//-- function
void load_database(FILE *fdb, database_t *DB)
{
  for (int i=0; i<NB_DPU_PER_RANK; i++) DB->db_prot[i] = create_db_prot(DB->nb_seq_dpu);
  int num_dpu = 0;
  protein_t *P;
  while ((P=fasta_get_protein(fdb))!=NULL)
    {
      int k = DB->db_prot[num_dpu]->nb_seq;
      DB->db_prot[num_dpu]->protein[k] = P;
      DB->db_prot[num_dpu]->nb_seq += 1;
      DB->db_prot[num_dpu]->nb_aa += P->length;
      num_dpu = (num_dpu + 1) % NB_DPU_PER_RANK;
    }
}

/*
//-- struct
typedef struct {
  uint64_t nb_aa;
  uint32_t nb_seq;
} db_info_t;
*/

//-- function
db_info_t *get_db_info(char *file_db_name)
{
  char *file_stat_name = (char *)malloc(sizeof(char)*(strlen(file_db_name)+10));
  int i;
  for (i=0; i<strlen(file_db_name); i++) file_stat_name[i] = file_db_name[i];
  file_stat_name[i] = '.';
  file_stat_name[i+1] = 'i';
  file_stat_name[i+2] = 'n';
  file_stat_name[i+3] = 'f';
  file_stat_name[i+4] = 'o';
  file_stat_name[i+5] = '\0';
  db_info_t *STAT = (db_info_t *) malloc(sizeof(db_info_t));
  if (access(file_stat_name, F_OK)==0)
    {
      FILE *ff = fopen(file_stat_name,"r");
      fread(STAT,sizeof(db_info_t),1,ff);
      fclose(ff);
    }
  else
    {
      FILE *ff = fopen(file_db_name,"r");
      protein_t *P;
      STAT->nb_aa = 0;
      STAT->nb_seq = 0;
      while ((P = fasta_get_protein(ff))!=NULL)
	{
	  STAT->nb_seq += 1;
	  STAT->nb_aa += P->length;
	  free_protein(P);
	}
      fclose(ff);
      ff = fopen(file_stat_name,"w");
      fwrite(STAT,sizeof(db_info_t),1,ff);
      fclose(ff);
    }
  free (file_stat_name);
  return STAT;
}

