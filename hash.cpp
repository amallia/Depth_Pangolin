/* hash.cpp */

#include <stdio.h>
#include <stdlib.h>
#include "hash.h"

/*********************************************************/
/* implementation of hash tables based on double hashing */
/*********************************************************/

 
/****************************************/
/* routines implementing the hash table */
/****************************************/

 
/****************************************/
/* routines implementing the hash table */
/****************************************/


/* To create hash table with at most limit elements.       */
/* Use exact=0 unless you know EXACTLY what you are doing. */
// void *initHash(int limit, int exact)
hashTable *initHash(int limit, int exact)

{
  hashTable *ht;
  int i;
  int primes [20] = { 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};

  ht = (hashTable *) malloc(sizeof(hashTable));
  if (exact == 1)
    ht->size = limit;
  else
  {
    for (i = 0; primes[i] < 1.35 * limit; i++);
    ht->size = primes[i];
  }
  ht->table = (short *) malloc(ht->size * sizeof(short));
  if (ht->table == NULL)  printf("Hash-internal error: Cannot alloc memory!\n");

  for (i = 0; i < ht->size; i++)  ht->table[i] = -1;
  return(ht);
}


/* To erase all elements in a hash table. */
/* Note: may leave unreferenced objects.  */
void clearHash(hashTable *ht)

{
  int i;

  for (i = 0; i < ht->size; i++)  ht->table[i] = -1;
}


/* To destroy a hash table.              */
/* Note: may leave unreferenced objects. */
void destroyHash(hashTable *ht)

{
  free(ht->table);
  free(ht);
}


/* To insert an element (pointer) into the hash table.   */
/* Returns 1 if key already exists, and 0 otherwise.     */
/* Note: if overwrite!=0 then existing item is replaced  */
/*   by new one; this may leave an unreferenced object.  */
int insertHash(hashTable *ht, int key, short elem, int overwrite, vector<int> & a)

{
  // printf("key to insert: %d ", key);
  // printf("elem to insert: %d \n", elem);
  int atry;
  unsigned int pos;

  // key = GETKEY(elem, a);
  atry = 0;
  pos = HASH(key, atry, ht->size);
  while ((ht->table[pos] != -1) && (GETKEY(ht->table[pos], a) != key))
  {
    atry++;
    if (atry > 100000)
    {
      printf("Hash-internal error: too many tries during insert!\n");
      return(0);
    }
    pos = HASH(key, atry, ht->size);
  }
    // work += atry+1;

  if (ht->table[pos] == -1)  
  {
    // printf("pos to insert: %d \n", pos);
    ht->table[pos] = elem;
    return(0);
  }
  else
  { 
    if (overwrite)  ht->table[pos] = elem;
    return pos+1;
    // return(1);
  }
}


/* to lookup an element in the hash table */
/* returns pointer to element if successful, and NIL otherwise */
short lookupHash(hashTable *ht, int key, vector<int>& a)
{
  int atry;
  unsigned int pos;

  atry = 0;
  pos = HASH(key, atry, ht->size);
  while ((ht->table[pos] != -1) && (GETKEY(ht->table[pos], a) != key))
  {
    atry++;
    if (atry > 100000)
    {
      printf("Hash-internal error: too many tries during lookup!\n");
      return(0);
    }
    pos = HASH(key, atry, ht->size);
  }

  printf("lookup, find at the pos: %d \n", pos);
  return(ht->table[pos]);
}