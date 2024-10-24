
#include "struct.h"


void write_align_format(FILE *falign, align_t *A, param_t *P)
{
  int i, j;
  int qs = A->query_start;
  int bs = A->db_start;
  uint8_t *Q = A->query_subseq;
  uint8_t *B = A->db_subseq;
  uint8_t *M = A->match;
  
  for (i=0; i<A->len; i+=60)
    {
      fprintf (falign,"Query %5d ",qs+i+1);
      for (j=i; j<i+60; j++)
        {
          if (j<A->len)
            {
              fprintf (falign,"%c",AA[Q[j]]);
            }
	  else
	    {
	      break;
	    }
        }
      fprintf (falign,"  %5d\n",qs+j);
      fprintf (falign,"            ");
      for (j=i; j<i+60; j++)
        {
          if (j<A->len)
	    {
	      fprintf (falign,"%c",M[j]);
	    }
        }
      fprintf (falign,"\n");
      fprintf (falign,"Sbjct %5d ",bs+i+1);
      for (j=i; j<i+60; j++)
        {
          if (j<A->len)
            {
              fprintf (falign,"%c",AA[B[j]]);
            }
	  else
	    {
	      break;
	    }
        }
      fprintf (falign,"  %5d\n\n",bs+j);
    }
}

void write_align(FILE *falign, param_t *P, query_t *Q)
{
  if (Q->nb_align == 0) return;
  align_t **AA = Q->align;
  protein_t *query = Q->protein;
  if (P->format == 6)
    {
      align_t *A;
      int nb_write_align = 0;
      for (int i=0; i<Q->nb_align; i++)
	{
	  A = AA[i];
	  fprintf (falign,"%s\t",A->query_com);
	  fprintf (falign,"%s\t",A->db_com);
	  double idt = (double) A->nb_ident; idt = (idt * 100) * (1.0); idt = idt/A->len;
	  fprintf (falign,"%2.2f\t",idt);
	  fprintf (falign,"%d\t",A->len);
	  fprintf (falign,"%d\t",A->len-A->nb_ident);
	  fprintf (falign,"0\t");
	  fprintf (falign,"%d\t",A->query_start+1);
	  fprintf (falign,"%d\t",A->query_start+A->len);
	  fprintf (falign,"%d\t",A->db_start+1);
	  fprintf (falign,"%d\t",A->db_start+A->len);
	  fprintf (falign,"%3.2e\t",A->evalue);
	  fprintf (falign,"%5.1f",A->bitscore);
	  fprintf (falign,"\n");
	  nb_write_align += 1;
	}
    }
  if (P->format == 0)
    {
      if (ftell(falign)==0)
        {
          fprintf (falign,"\nDatabase: %s\n", P->db_file_name);
        }
      align_t *A = AA[0];
      fprintf (falign,"Query= %s\n\nLength=%d\n\n",A->query_com,query->length);
      fprintf (falign,"                                                                       Raw     Bit\n");
      fprintf (falign,"Sequences producing significant alignments:                            Score   Score    E-Value\n\n");
      int nb_write_align = 0;
      for (int i=0; i<Q->nb_align; i++)
	{
	  A = AA[i];
	  fprintf (falign," ");
	  int i = 0;
	  while ((A->db_com[i]!='\0')&&(i<67)) { fprintf (falign,"%c",A->db_com[i]); i += 1; }
	  fprintf (falign," ");
	  while (i<68) { fprintf (falign,"."); i += 1; }
	  fprintf (falign," %5d   %5.1f   %3.2e\n",A->score,A->bitscore,A->evalue);
	  nb_write_align += 1;
	}
      fprintf (falign,"\n\n");
      nb_write_align = 0;
      for (int i=0; i<Q->nb_align; i++)
	{
	  A = AA[i];
	  fprintf (falign,"%s\n\n",A->db_com);
	  fprintf (falign,"            Raw Score = %5d, Bit Score = %5.2f, E-Value = %3.2e\n",A->score,A->bitscore,A->evalue);
	  fprintf (falign,"            Identities = %d/%d (%d%c), ",A->nb_ident,A->len,(A->nb_ident*100)/A->len,37);
	  fprintf (falign,"Positives = %d/%d (%d%c), ",A->nb_pos,A->len,(A->nb_pos*100)/A->len,37);
	  fprintf (falign,"Gaps = %d/%d (%d%c)\n\n",0,A->len,0,37);
	  write_align_format(falign,A,P);
	  nb_write_align += 1;
	}
    }
}

//-- function
void save_aligns(syncbuf_t *SBAlign, param_t *P, FILE *falign)
{
  query_t *Q;
  int nbq = 1;
  while ((Q = syncbuf_get(SBAlign))!=NULL)
    {
      write_align(falign,P,Q);
      free_query(Q);
      if (P->verbose) { printf ("  - %d%c",nbq,13); fflush(stdout); }
      nbq += 1;
    }
  if (P->verbose) printf ("\n");
}
