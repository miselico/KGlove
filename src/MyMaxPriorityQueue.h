
#ifndef MY_MAXPRIORITYQUEUE_H_
#define MY_MAXPRIORITYQUEUE_H_

#include <vector>
#include <unordered_map>

/**
 * Adapted from Snap/glib priorityqueue
 */

template<class VAL>
class MyMaxPriorityQueue {

private:
	std::vector<double> Priorities;
	std::unordered_map<VAL, unsigned int> ValToIndex;
	std::vector<VAL> IndexToVal;

	unsigned int Parent(unsigned int i);
	unsigned int Left(unsigned int i);
	unsigned int Right(unsigned int i);

	void Swap(unsigned int i, unsigned int j);

	// If the max-heap invariant is satisfied except for index i possibly being smaller than a child, restore the invariant.
	void MaxHeapify(unsigned int i);
public:
	void Insert(const VAL& X, double Priority);

	void SetPriority(const VAL& X, double NewPriority);

	double GetPriority(const VAL& X);

	float GetMaxPriority();

	int PopMax();

	int Size();

	bool IsEmpty();
};

template<class VAL> unsigned int MyMaxPriorityQueue<VAL>::Parent(unsigned int i) {
	return (i + 1) / 2 - 1;
}

template<class VAL> unsigned int MyMaxPriorityQueue<VAL>::Left(unsigned int i) {
	return i * 2 + 1;
}

template<class VAL> unsigned int MyMaxPriorityQueue<VAL>::Right(unsigned int i) {
	return i * 2 + 2;
}

template<class VAL> void MyMaxPriorityQueue<VAL>::Swap(unsigned int i, unsigned int j) {
	std::swap(Priorities[i], Priorities[j]);
	std::swap(IndexToVal[i], IndexToVal[j]);

	ValToIndex[IndexToVal[i]] = i;
	ValToIndex[IndexToVal[j]] = j;
}

// If the max-heap invariant is satisfied except for index i possibly being smaller than a child, restore the invariant.
template<class VAL> void MyMaxPriorityQueue<VAL>::MaxHeapify(unsigned int i) {
	unsigned int largest = i;
	if (Left(i) < Priorities.size() && Priorities[Left(i)] > Priorities[largest]) {
		largest = Left(i);
	}
	if (Right(i) < Priorities.size() && Priorities[Right(i)] > Priorities[largest]) {
		largest = Right(i);
	}
	if (largest != i) {
		Swap(i, largest);
		MaxHeapify(largest);
	}
}

template<class VAL> void MyMaxPriorityQueue<VAL>::Insert(const VAL& X, double Priority) {
	ValToIndex[X] = IndexToVal.size();
	IndexToVal.push_back(X);
	// Priorities.Add(-INFINITY);
	Priorities.push_back(-1.0);
	SetPriority(X, Priority);
}

template<class VAL> void MyMaxPriorityQueue<VAL>::SetPriority(const VAL& X, double NewPriority) {
	auto i_iter = ValToIndex.find(X);
	if (i_iter == ValToIndex.end()) {
		Insert(X, NewPriority);
	} else {
		int i = i_iter->second;
		if (NewPriority >= GetPriority(X)) {
			Priorities[i] = NewPriority;
			while (i > 0 && Priorities[i] > Priorities[Parent(i)]) {
				Swap(i, Parent(i));
				i = Parent(i);
			}
		} else {
			Priorities[i] = NewPriority;
			MaxHeapify(i);
		}
	}
}

template<class VAL> double MyMaxPriorityQueue<VAL>::GetPriority(const VAL& X) {
	auto index = ValToIndex.find(X);
	if (index != ValToIndex.end()) {
		return Priorities[index->second];
	} else {
		return 0.0;
	}
}

template<class VAL> float MyMaxPriorityQueue<VAL>::GetMaxPriority() {
	assert(Size() > 0 && "Attempt to query max priority of empty priority queue.");
	return Priorities[0];
}

template<class VAL> int MyMaxPriorityQueue<VAL>::PopMax() {
	assert(Size() > 0 && "Attempt to query max priority of empty priority queue.");
	int maxVal = IndexToVal[0];
	Swap(0, Priorities.size() - 1);
	Priorities.pop_back();
	IndexToVal.pop_back();
	ValToIndex.erase(maxVal);

	MaxHeapify(0);
	return maxVal;
}

template<class VAL> int MyMaxPriorityQueue<VAL>::Size() {
	return Priorities.size();
}

template<class VAL> bool MyMaxPriorityQueue<VAL>::IsEmpty() {
	return Size() == 0;
}

#endif /*MY_MAXPRIORITYQUEUE_H_*/
