
#include "struct.h"

/*
//-- struct
typedef struct {
  sem_t mutex;
  sem_t empty;
  sem_t full;
  char *name;
  int nb_cons;
  int nb_end;
  int nb_elt;;
  int nb_empty;
  int nb_full;
  int size_buf;
  int index_write;
  int index_read;
  query_t **buf;
} syncbuf_t;
*/

//-- function
syncbuf_t *create_syncbuf(char *name, int size_buffer, int nb_consumer)
{
  syncbuf_t *SB = (syncbuf_t *) malloc(sizeof(syncbuf_t));
  sem_init(&SB->mutex,0,1);
  sem_init(&SB->empty,0,size_buffer);
  sem_init(&SB->full,0,0);
  SB->name = name;
  SB->nb_cons = nb_consumer;
  SB->nb_end = 0;
  SB->index_write = 0;
  SB->index_read = 0;
  SB->nb_elt = 0;
  SB->nb_empty = 0;
  SB->nb_full = 0;
  SB->size_buf = size_buffer;
  SB->buf = (query_t **) malloc(size_buffer*sizeof(query_t *));
  return SB;
}

//-- function
void syncbuf_free (syncbuf_t *SB)
{
  free (SB->buf);
  free (SB);
}

//-- function
void syncbuf_insert(syncbuf_t *SB, query_t *Q)
{
  sem_wait (&SB->empty);
  sem_wait (&SB->mutex);
  SB->buf[SB->index_write] = Q;
  SB->index_write = (SB->index_write + 1) % SB->size_buf;
  SB->nb_elt += 1;
  if (SB->nb_elt == SB->size_buf) SB->nb_full += 1;
  sem_post (&SB->mutex);
  sem_post (&SB->full);
}

//-- function
void *syncbuf_get(syncbuf_t *SB)
{
  query_t *Q;
  sem_wait (&SB->full);
  sem_wait (&SB->mutex);
  Q = SB->buf[SB->index_read];
  SB->index_read = (SB->index_read + 1) % SB->size_buf;
  SB->nb_elt -= 1;
  if (SB->nb_elt == 0) SB->nb_empty += 1;
  sem_post (&SB->mutex);
  sem_post (&SB->empty);
  return Q;
}

//-- function
int syncbuf_end(syncbuf_t *SB)
{
  SB->nb_end += 1;
  if (SB->nb_end < SB->nb_cons)
    {
      syncbuf_insert(SB,NULL);
      return 0;
    }
  else
    {
      return 1;
    }
}

//-- function
void syncbuf_printstat(syncbuf_t *SB)
{
  printf ("  - buffer %s size=%2d, #empty=%4d, #full=%4d\n",SB->name, SB->size_buf, SB->nb_empty, SB->nb_full);
}
