/************************************************/
/* definitions needed for hash function package */
/************************************************/

#define H1  2197
#define H2  1331
#define HASH(k,i,m) (unsigned int)((( ((long long int) k) * H1 +  \
   ((long long int) i) * (1+(((long long int) k)*H2)%((long long int)(m-1))))\
                          ) % (long long int)(m))
/*
#define HASH(k,i,m) (unsigned int)(( ((long long int) k) * H1 +  \
                       ((long long int) i) * ((1+((long long int) k)*H2)%((long long int)(m-1))))\
                          ) % (long long int)(m)
*/

#include <vector>

using namespace std;


typedef struct {
  int size;
  short * table; 
} hashTable;


// int work = 0; 



/* define here how to get integer key from an element         */
/* NOTE: needs to be adapted to the structure of the elements */
// #define GETKEY(e) (*((int *)(e)))



// hashTable *initHash(int limit, int exact);

hashTable *initHash(int limit, int exact);

void clearHash(hashTable *ht);

void destroyHash(hashTable *ht);

int insertHash(hashTable *ht, int key, short elem, int overwrite, vector<int> & a);

short lookupHash(hashTable *ht, int key, vector<int> & a);

template<typename T>
int GETKEY(short e, vector<T> & a){
  // return *(a+e);
	if(e<=a.size()){
  		return a.at(e);
	}else{
		return -1;
	}
}