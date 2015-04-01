#ifndef HEAP_H
#define HEAP_H

#include <vector>

using namespace std;

//
// InversionHeap is designed to obtain the inversion with lowest cost quickly
// Two adjacent variables form an inversion
// if their relative order is different from that in the final permutation
//

typedef struct Inversion
{
	// An inversion is formed by the variables at level and level+1
	int level;
	long cost;
  Inversion(int l, long c) : level(l), cost(c) {}
} Inversion;

class InversionHeap
{
private:
	static const int NOT_IN_HEAP=-1;

	static inline size_t left(size_t i)
	{
		return i*2+1;
	}

	static inline size_t right(size_t i)
	{
		return i*2+2;
	}

	static inline size_t parent(size_t i)
	{
		return (i-1)/2;
	}

	static inline bool lessthan(Inversion& left, Inversion& right)
	{
		return left.cost<right.cost;
	}

	vector<Inversion> _heap;
	vector<int> _indices;

	bool inline has_left(size_t index)
	{
		return left(index)<_heap.size();
	}

	bool inline has_right(size_t index)
	{
		return right(index)<_heap.size();
	}

	void percolate_down(int level)
	{
		assert(is_in_heap(level));
		int i=_indices[level];
		Inversion inversion=_heap[i];
		while(has_left(i)){
			int min_index = left(i);
			if (has_right(i)){
				if (lessthan(_heap[right(i)], _heap[min_index])) {
					min_index = right(i);
				}
			}
			if (lessthan(_heap[min_index], inversion)) {
				_indices[_heap[min_index].level]=i;
				_heap[i]=_heap[min_index];
				i=min_index;
			}
			else{
				break;
			}
		}
		_indices[level]=i;
		_heap[i]=inversion;
	}

	void percolate_up(int level)
	{
		assert(is_in_heap(level));
		int i=_indices[level];
		Inversion inversion=_heap[i];
		int p=parent(i);
		while (i>0 && lessthan(inversion, _heap[p])){
			_indices[_heap[p].level]=i;
			_heap[i]=_heap[p];
			i=p;
			p=parent(i);
		}
		_indices[level]=i;
		_heap[i]=inversion;
	}

public:
	InversionHeap(int var_size)
	{
		_indices.assign(var_size+1, NOT_IN_HEAP);
		_heap.reserve(var_size);
	}

	// If the level is in the heap, update the cost
	void push(int level, long cost)
	{
		if(!is_in_heap(level)){
			_indices[level]=_heap.size();
			_heap.push_back(Inversion(level, cost));
			percolate_up(level);
		}
		else {
			long old_cost=_heap[_indices[level]].cost;
			_heap[_indices[level]].cost=cost;
			if(old_cost<cost) {
				percolate_down(level);
			}
			else if(old_cost>cost) {
				percolate_up(level);
			}
		}
	}

	void pop()
	{
		int top=_heap[0].level;
		_heap[0]=_heap.back();
		_indices[_heap[0].level]=0;
		_indices[top]=NOT_IN_HEAP;
		_heap.pop_back();
		if(_heap.size()>1) {
			percolate_down(_heap[0].level);
		}
	}

	int top()
	{
		return _heap[0].level;
	}

	bool is_in_heap(int level) const
	{
		return _indices[level]!=NOT_IN_HEAP;
	}

	bool empty() const
	{
		return _heap.empty();
	}
};


#endif
