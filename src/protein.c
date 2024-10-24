
#include "struct.h"

/*
//-- struct
typedef struct {
  char *com;
  uint8_t *seq;
  uint32_t length;
} protein_t;
*/

//-- function
void free_protein(protein_t *P)
{
  if (P==NULL) return;
  if (P->com != NULL) free (P->com);
  if (P->seq != NULL) free (P->seq);
  free (P);
}

//-- extern
char AA[32] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*','-','#','$','@','!','+','&','%','?','!','='}; 

//-- extern
int CODE_AA [256] = 
  { 
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /*   0 -  15 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /*  16 -  31 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /*  32 -  47 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /*  48 -  63 */
    20,  0, 20,  4,  3,  6, 13,  7,  8,  9, 20, 11, 10, 12,  2, 20,           /*  64 -  79 */
    14,  5,  1, 15, 16, 20, 19, 17, 20, 18, 20, 20, 20, 20, 20, 20,           /*  80 -  95 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /*  96 - 111 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,           /* 112 - 127 */
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
  };

//-- function
uint8_t *filter_complexity(protein_t *P)
{
  int WS = 20;
  int X = 0;
  uint8_t *cplx = (uint8_t *) calloc(P->length,sizeof(uint8_t));
  if (P->length > WS+5)
    {
      int *T = (int *) calloc(32, sizeof(int));
      int8_t aa;
      for (int i=0; i<WS; i++)
	{
	  aa = P->seq[i];
	  T[aa] += 1;
	  if (T[aa]==1) X += 1;
	}
      for (int i=WS; i<P->length; i++)
	{
	  aa = P->seq[i-WS];
	  T[aa] -= 1;
	  if (T[aa]==0) X -= 1;
	  aa = P->seq[i];
	  T[aa] += 1;
	  if (T[aa]==1) X += 1;
	  if (X<5)
	    for (int j=i-WS; j<i; j++) cplx[j] = 1;
	}
      free(T);
    }
  return cplx;
}



// read sequence from FASTA format
// '>' is not kept inside the comments

//-- function
protein_t *fasta_get_protein(FILE *fastafile)
{
  if (feof(fastafile)) return NULL;
  protein_t *P = (protein_t *)malloc(sizeof(protein_t));
  // read comment
  uint64_t start_com = ftell(fastafile);
  char c;
  int l = 0;
  while ((c=fgetc(fastafile))!='\n') l += 1;
  P->com = (char *) malloc(sizeof(char)*(l));
  fseek(fastafile,start_com,SEEK_SET);
  c=fgetc(fastafile);
  l = 0;
  while ((c=fgetc(fastafile))!='\n') { P->com[l] = c; l += 1; }
  P->com[l] = '\0';
  uint64_t start_seq = ftell(fastafile);
  l = 0;
  while ((c=fgetc(fastafile))!='>')
    {
      if (c==EOF) break;
      if (c!='\n') l += 1;
    }
  P->seq = (uint8_t *) malloc(sizeof(uint8_t)*(l+1));
  fseek(fastafile,start_seq,SEEK_SET);
  l = 0;
  while ((c=fgetc(fastafile))!='>')
    {
      if (c==EOF) break;
      if (c!='\n') { P->seq[l] = CODE_AA[(int)c]; l += 1; }
    }
  P->seq[l] = '\0';
  P->length = l;
  if (c!=EOF) fseek(fastafile,ftell(fastafile)-1,SEEK_SET);
  return P;
}


//-- function
void display_protein(protein_t *P)
{
  printf (">%s\n",P->com);
  int k=1;
  for (int i=0; i < P->length; i++)
    {
      printf ("%c",AA[P->seq[i]]);
      if (k%60 == 0) printf ("\n");
      k += 1;
    }
  if ((k-1)%60 != 0) printf ("\n");
}
