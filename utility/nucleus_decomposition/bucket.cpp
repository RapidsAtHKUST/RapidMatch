#include "bucket.h"

/******* ARRAY BUCKET ********/

Naive_Bucket_element::Naive_Bucket_element()
: next(NULL), prev(NULL)
{}

Naive_Bucket::Naive_Bucket()
: max_value(0), current_min_value(1), elements(NULL), buckets(NULL), values(NULL), nb_elements(0)
{}

Naive_Bucket::~Naive_Bucket() {
	if (buckets != NULL)
		free (buckets);
	if (elements != NULL)
		free (elements);
	if (values != NULL)
		free (values);
}

void Naive_Bucket::Free () {
	free (buckets);
	buckets = NULL;
	free (elements);
	elements = NULL;
	free (values);
	values = NULL;
}

void Naive_Bucket::Initialize(int max_v, int nb_element) {
	int i;
	max_value = max_v;
	buckets = (Naive_Bucket_element **) malloc(sizeof(Naive_Bucket_element *) * (max_value+1));
	elements = (Naive_Bucket_element *) malloc(sizeof(Naive_Bucket_element) * nb_element);
	values = (int *) malloc(sizeof(int) * nb_element);
	nb_elements = nb_element;
	if (buckets == NULL || elements == NULL || values == NULL) {
		free(values);
		free(buckets);
		free(elements);
	}
	else {
		for (i = 0; i <= max_value; i++)
			buckets[i] = NULL;
		for (i = 0; i < nb_element; i++) {
			elements[i].prev = NULL;
			elements[i].next = NULL;
		}
	}
	current_min_value = max_value + 1;
}

void Naive_Bucket::Insert (int id, int value) {
	values[id] = value;
	elements[id].prev = NULL;
	elements[id].next = buckets[value];
	if (buckets[value] != NULL)
		buckets[value]->prev = &(elements[id]);
	else if (current_min_value > value)
		current_min_value = value;
	buckets[value] = &(elements[id]);
}

int Naive_Bucket::PopMin(int* id, int* ret_value) {
	for (; current_min_value <= max_value; current_min_value++) {
		if (buckets[current_min_value] != NULL) {
			*id = buckets[current_min_value] - elements; // pointer arithmetic. finds the index of element that buckets[current_min_value] points to
			buckets[current_min_value] = buckets[current_min_value]->next; // adjust the pointer to the new head of the list that has same degree elements
			if (buckets[current_min_value] != NULL)
				buckets[current_min_value]->prev = NULL;
			*ret_value = values[*id];
			values[*id] = -1;
			return 0;
		}
	}
	return -1; // if the bucket is empty
}

int Naive_Bucket::CurrentValue(int id) {
	return values[id];
}

void Naive_Bucket::DecVal(int id) {
	int old_value = values[id];
	// adjust the prev and next pointers
	if (elements[id].prev == NULL)
		buckets[old_value] = elements[id].next;
	else
		elements[id].prev->next = elements[id].next;
	if (elements[id].next != NULL)
		elements[id].next->prev = elements[id].prev;
	Naive_Bucket::Insert(id, values[id]-1);
	return;
}

