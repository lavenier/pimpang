
#include "struct.h"


/*
//-- struct
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
*/

int GetCPUCount()
{
 cpu_set_t cs;
 long nproc = sysconf(_SC_NPROCESSORS_ONLN);
 CPU_ZERO(&cs);
 sched_getaffinity(0, sizeof(cs), &cs);

 int count = 0;
 for (int i = 0; i < nproc; i++)
 {
  if (CPU_ISSET(i, &cs))
   count++;
 }
 return count;
}

//-- function
param_t *get_parameter (int argc, char *argv[])
{
  param_t *PARAM = (param_t *) malloc(sizeof(param_t));
  
  PARAM->pim_mode = 0;
  PARAM->evalue = DEFAULT_EVALUE;
  PARAM->db_file_name = "none";
  PARAM->query_file_name = "none";
  PARAM->align_file_name = "none";
  PARAM->submat_name = DEFAULT_SUB_MATRIX;
  PARAM->lambda = DEFAULT_LAMBDA;
  PARAM->K = DEFAULT_K;
  PARAM->format = 6;
  PARAM->verbose = 0;
  PARAM->nb_threads = GetCPUCount();
  PARAM->nb_ranks = 0;

  char c;
  int help = 0;
  while ((c = getopt(argc, argv, "p:d:q:o:e:f:m:vh")) >= 0)
    {
           if (c == 'm') { PARAM->submat_name      = optarg; }
      else if (c == 'd') { PARAM->db_file_name     = optarg; }
      else if (c == 'q') { PARAM->query_file_name  = optarg; }
      else if (c == 'o') { PARAM->align_file_name  = optarg; }
      else if (c == 'e') { PARAM->evalue           = atof(optarg); }
      else if (c == 'p') { PARAM->nb_ranks         = atoi(optarg);  PARAM->pim_mode = 1; }
      else if (c == 'f') { PARAM->format           = atoi(optarg); }
      else if (c == 'v') { PARAM->verbose          = 1; }
      else if (c == 'h') { help                    = 1; }
      else exit(0);
    }
  if (help)
    {
      printf ("USAGE : %s [options] -d db -q query -o results\n",argv[0]);
      printf ("options\n");
      printf ("-m <string> : substitution matrix name [blosum90-80-62-50-45 pam30-70-250](default = blosum62)\n");
      printf ("-d <string> : data base file name [mandatory]\n");
      printf ("-q <string> : query file name [mandatory]\n");
      printf ("-o <string> : alignment file name [mandatory]\n");
      printf ("-e <float>  : evalue (default = %3.2e)\n",DEFAULT_EVALUE);
      printf ("-p <int>    : activate PIM mode; number of DPU ranks (max = %d)\n",MAX_NB_RANKS);
      printf ("-f <int>    : format = 0 | 6 (default = 6)\n");
      printf ("-v          : verbose mode\n");
      printf ("-h          : help\n");
      printf ("\n");
      exit (1);
    }

  if (PARAM->nb_ranks > MAX_NB_RANKS)
    {
      PARAM->nb_ranks = MAX_NB_RANKS;
    }
  if (strcmp(PARAM->db_file_name,"none")==0)
    {
      fprintf (stderr,"\nERROR: db file name not provided (-d)\nexit...\n\n"); exit (0);
    }
  if (strcmp(PARAM->query_file_name,"none")==0)
    {
      fprintf (stderr,"\nERROR: query file name not provided (-q)\nexit...\n\n"); exit (0);
    }
  if (strcmp(PARAM->align_file_name,"none")==0)
    {
      fprintf (stderr,"\nERROR: alignment file name not provided (-o)\nexit...\n\n"); exit (0);
    }
  return PARAM;
}

