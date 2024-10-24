
#include "struct.h"

// transform raw score to bit score
//-- function
double bitscore(param_t *P, int32_t score)
{
  double bs;
  double sc = (double) score;
  bs = P->lambda * sc;
  bs = bs - logf(P->K);
  double two = 2.0;
  bs = bs / log(two);
  return bs;
}

// compute evalue from a score
//-- function
double score2evalue(param_t *P, int query_len, uint64_t nb_aa_db, int score)
{
  double QL = (double) query_len;
  double NBAA = (double) nb_aa_db;
  return NBAA * QL * P->K * exp((-P->lambda * score)); 
}

// compute score from e-value
//-- function
uint32_t evalue2score(param_t *P, int query_len, uint64_t nb_aa_db)
{
  double QL = (double) query_len;
  double NBAA = (double) nb_aa_db;
  double score = (-(log(P->evalue/(NBAA * QL * P->K))/P->lambda));
  return (uint32_t) score;
}
