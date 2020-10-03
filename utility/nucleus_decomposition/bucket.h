#ifndef __BUCKET__H
#define __BUCKET__H

#include <vector>
#include <unordered_map>
#include <stdlib.h>

struct Naive_Bucket_element
{
	Naive_Bucket_element * prev;
	Naive_Bucket_element * next;
	Naive_Bucket_element();
};

struct Naive_Bucket
{
	Naive_Bucket_element **buckets; /* actual pointers to bucket heads */
	Naive_Bucket_element *elements; /* for direct access to bucket elements elements[id] is the id-th element */
	int nb_elements;
	int max_value;
	int *values; /* needed for update, in case bucket head changed. */
	int current_min_value;

	Naive_Bucket();
	~Naive_Bucket();
	void Initialize(int max_value, int nb_element);
	void Free (); /* value == INT_MAX means not present in bucket */
	void Insert(int id, int value);
	void Update(int id, int new_value);
	void DecVal(int id);
	int PopMin(int* ret_id, int* ret_value); /* returns -1 if empty */
	int CurrentValue(int id);
};

#endif
