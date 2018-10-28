###
## Comparative performance evaluation of sorting algorithms

Author: Kejie Zhang
November 28th,2017



## Abstract

This article presents a comprehensive performance evaluation and comparison of four kinds of sorting algorithms: Insertion Sort, Merge Sort, Quick Sort and Optimized Quick Sort. Using different sizes and sequential circumstances to calculate the time cost and number of comparisons in order to evaluate and compare the three basic and one optimized sorting algorithms. This article considers two metrics: the number of comparisons and the time cost, and four different initial conditions in measuring the efficiency of difference algorithms.

## Introduction

Sorting problem is one of the fundamental aspects of algorithms. We take an array as input, perform a sequence of specified operations on it and output a sorted array. There are different kinds of sorting algorithms, we are focusing on the specific ones: Comparison Sorts. A Comparison Sort is a sorting algorithm where the final order is determined only by comparisons between the input elements. We compare elements at each step of the algorithm to determine whether one element should be on the right or left of another element. Using comparison sort, we can quickly and easily to solve out the most unsorted integer arrays.

Comparison sorts are straightforward and obvious: we compare two elements step by step and generate a sequence of ascendant or descendant elements. There are many ways to choose which two elements to compare each other, thus there are different methods to implement each algorithms. Some basic comparison sorting algorithms are like Selection Sort, Bubble Sort, Insertion Sort, Merge Sort, Heap Sort, Quick Sort, Radix Sort, Shell Sort and etc. In this article, we only consider four algorithms: Insertion Sort, Merge Sort, Quick Sort and Optimized Quick Sort. We will input four different initial conditions of arrays in different sizes. We will measure time complexity in worse-case, average-case and best-case by drawing 2-D graphs.

Comparison sorts can be viewed abstractly in terms of decision trees. A decision tree is a full binary tree that represents the comparisons between elements that are performed by a particular sorting algorithm operating on an input of a given size. The execution of the sorting algorithm corresponds to tracing a path from the root of the decision tree to a leaf. From the decision tree definition, we can know that the lower bound of all comparison sorts is O(nlogn). The lower bound means to use the best algorithm to run in the worst-case circumstance. So comparison sorts cannot run faster than O(nlogn).

Although comparison sorts all focus in one goal: to generate a sorted list of elements, however, the different way they implement will result different time complexity. We will use the Big O to discuss their time complexity. And sometimes the input array is in different initial situation: a random initial order, a sorted order, a reversed order and an almost the same with few unique keys order. We can see that different algorithms performance different efficiency in these four circumstances. Actually we can have a conclusion in the later that there is no best sorting algorithm, only the best option for the specific circumstance.

## Stimulation Experiment

The stimulation environment is run in a 2.7GHz Intel Core i7 running MacOS High Sierra. The IDE is Xcode 9.1 and the algorithms are written in C++. Five different array sizes are used: 10000, 20000, 30000, 40000 and 50000. I choose these sizes of numbers because they can show us a direct All the arrays are generated in random. All the point in figures are statistically meaningful, with 50 times run first.

**Insertion sort**

Insertion sort is definitely not the most efficient algorithm but it is appraised by its simplicity. It builds a final sorted array one element at a time. It iterates through an input array and removes one element per iteration, finds the place the element belongs in the array, and then places it there. Since it is very easy to implement and adequately efficient for small number of elements, it is useful for people to use in certain kinds of circumstances. Usually, when the number of Array is less than 100, we can use insertion sort to complete the sorting issue. Besides, insertion sort is powerful when the sequence is almost sorted. In best-case, the time complexity of it will beO(n). However, the running time of insertion sort is O(n2), in average-case and worst-case. So if we want to sort a random unsorted array with a huge size, the performance of insertion sort will become extremely bad.

| ComparisonNum | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 24991700 | 99993800 | 225016000 | 400154000 | 624450000 |
| Sorted | 23636.4 | 47290.6 | 70899.8 | 94537.2 | 118154 |
| Reversed | 49986800 | 199974000 | 449960000 | 799947000 | 1249930000 |
| Duplicated | 24857200 | 99731400 | 224445000 | 399500000 | 624691000 |

| Time | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 0.125299 | 0.533186 | 1.20039 | 2.02659 | 3.16744 |
| Sorted | 0.00032336 | 0.00061656 | 0.00091846 | 0.00117968 | 0.00147514 |
| Reversed | 0.262502 | 1.04695 | 2.4269 | 4.28468 | 6.6767 |
| Duplicated | 0.125078 | 0.50346 | 1.2092 | 2.08798 | 3.38729 |

**Merge Sort**

Merge sort focuses on how to merge two pre-sorted arrays together such that the resulting array is also sorted. Unlike insertion sort, merge sort has a time complexity of O(nlog(n)) in best-case, average-case, and worst-case. So it is efficient in all kinds of initial situations. The limit of this algorithm is its space inefficiency. Merge sort actually save its time by wasting more space, it is a kind of time-space trade off. A lot of temporary arrays have to be created and many copying of elements is involved during the operation, so the worst-case space complexity is O(n). Despite of this, the merge sort is still very proficient: it is stable and it is a truly O(nlog(n)) algorithm.

| ComparisonNum | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 120456 | 260894 | 408605 | 561798 | 718195 |
| Sorted | 73927.1 | 157856 | 244565 | 335767 | 426356 |
| Reversed | 73258.5 | 156498 | 243473 | 332996 | 425308 |
| Duplicated | 120388 | 260837 | 408505 | 561736 | 718113 |

| Time | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 0.00221032 | 0.00448748 | 0.00673522 | 0.0098256 | 0.0121162 |
| Sorted | 0.00135708 | 0.00272294 | 0.00414936 | 0.00575542 | 0.00715248 |
| Reversed | 0.00136524 | 0.00286246 | 0.004544 | 0.00574088 | 0.0075217 |
| Duplicated | 0.00192388 | 0.00420078 | 0.00670894 | 0.00862678 | 0.0112535 |

**Quick Sort**

Quick sort is one of the most widely used sorting algorithms. Indeed, it has plenty virtue: it is in-place, it is really fast since its inner loop is very short, and its typical time complexity is O(nlog(n)). However, quick sort has a fatal problem, it will become very slow when it is in worst-case situation, having a O(n2) time complexity. But the worst-case is not common, if the array contains random elements without having too many duplicates, using quick sort is still a wise choice.

| ComparisonNum | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 77722.4 | 170898 | 267404 | 366227 | 471595 |
| Sorted | 22041800 | 88270700 | 198606000 | 352671000 | 550918000 |
| Reversed | 8845250 | 35380700 | 79795200 | 141783000 | 221526000 |
| Duplicated | 297604 | 609006 | 932011 | 1247840 | 1565780 |

| Time | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 0.00185642 | 0.0039663 | 0.00575812 | 0.0080376 | 0.0101239 |
| Sorted | 0.211133 | 0.845443 | 1.92674 | 3.45148 | 5.19816 |
| Reversed | 0.0957615 | 0.388007 | 0.916991 | 1.58424 | 2.48285 |
| Duplicated | 0.0036941 | 0.00724664 | 0.0112629 | 0.0154369 | 0.0192974 |

**Optimized Quick Sort**

Considering quick sort as a very efficient and popular algorithm, there are some optimizations to improve its performance. In this article, there are two different kinds of quick sort. One is a plain original one, implemented as the most basic method, the other is an optimized version. The optimized quick sort adopts some optimizations: picking the median as pivot, using two-way partitioning, using insertion sort when the sub array is in small size, and using tail recursion.

1. 1)Picking median as pivot:

The common quick sort is used to choose the rightest or leftist element of the array as a pivot. When it chooses a biggest or smallest one as pivot, it will make the quick sort become insertion sort and thus reduce its performance dramatically. One solution is to choose the pivot randomly, so the probability to choose the extremely big or small element will reduce to lowest. Indeed, it is useful, however, there is still a chance to let the pivot be the worst. So we can compare three elements-right, middle, and left ones, and choose the median one as the pivot. So the algorithm will never choose the biggest or smallest element.

1. 2)Two-way partitioning

If the array has many duplicate elements, we can find that basic quick sort will become nearly O(n2), because when we are in the partition process, we put all the elements that are equal to pivot on one side of the array, so in those special initial situations the performance will become bad. What can we do to optimize it? The solution is to use two-way partitioning, that is, to compare the pivot both from two sides of the array instead of one side. This will put the duplicated elements on both side of the array uniformly.

1. 3)Using insertion sort when the sub array is in small size

With the partitions are going, the sub array&#39;s size go down and we can use insertion sort in replace. We can define a number, so that when the size of sub an array is less than it, we will change to sort them in insertion sort.

1. 4)Using tail recursion

Quick sort needs extra space to recursive function calls, usually in worst case it needsO(n) space. Because that the recursive function appears after partition, so we can optimize it by using tail elimination. We convert it to recall only one recursive and reduce the extra space into O(logn) in worst-case.

| ComparisonNum | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 111618 | 236231 | 366003 | 502643 | 639969 |
| Sorted | 158953 | 347155 | 546121 | 747519 | 966230 |
| Reversed | 148572 | 325623 | 516357 | 705766 | 911619 |
| Duplicated | 64567.5 | 143594 | 229149 | 315132 | 404745 |

| Time | 10000 | 20000 | 30000 | 40000 | 50000 |
| --- | --- | --- | --- | --- | --- |
| Random | 0.00203838 | 0.00397546 | 0.00653634 | 0.00826116 | 0.0111022 |
| Sorted | 0.00128728 | 0.00257056 | 0.00382754 | 0.00514782 | 0.00628898 |
| Reversed | 0.00130454 | 0.00250766 | 0.00380722 | 0.0049461 | 0.00642398 |
| Duplicated | 0.00166144 | 0.00332054 | 0.00516452 | 0.00704436 | 0.00889354 |

#### Conclusion

From the previous result we can get the conclusions as follows:

1.
1.For the random array, the insertion sort performances the worst. Since it is an O(n2) algorithm, it will cost much more time and comparison times to sort. While the merge sort, quick sort and optimized quick sort performance pretty well, we can see that they have a O(nlogn) time-complexity.
2.
2.For the sorted array, it is the worst case for quick sort and we can see how bad it is, it becomes an O(n2) algorithm and it is really slow. Insertion sort, merge sort and optimized quick sort are still O(nlogn) algorithm. On the contrary, the insertion sort behaves pretty well in sorted situation. At the mean time, we can see that the solutions we made to optimized quick sort is obvious and efficient, the number of comparisons and time cost reduce dramatically.
3.
3.For the reversed sort, it is the worst case for quick sort and we can see it has the same time complexity with insertion sort, both are O(n2). They behave very bad in this situation. In other way, merge sort and optimized quick sort performance very well, they have the same complexity of O(nlogn).
4.
4.For the duplicated array, the insertion sort and quick sort performs bad, they are both O(n2) algorithm. And merge sort and the optimized quick sort still behave as O(nlogn).

We can see that insertion sort is good for sorted array, however, since this situation is so rare in practice, we do not want to use it in big size of input array. Merge sort is stable and trustful, it will remain in O(nlogn) in best, average and worst case. Quick sort is not stable and trustful, it can perform well in average case, however, will performs really bad in worst-case, which is sorted, reversed and numerous duplications.

The efficiency of optimized quick sort is very tremendous, it remains O(nlogn) in random, sorted and reversed initial situations. It reduces time in numerous duplication situation, but it still remains O(n2), but it is quite practice and enough to deal with most cases.