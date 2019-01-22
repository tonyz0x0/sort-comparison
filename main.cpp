/*********************************
 * Algorithm-Comparison
 * 
 * DESCRIPTION: A collection of various sorting algorithms and quantity evaluation methods
 *
 * Author: Kejie Zhang
 * LAST UPDATED: 01/20/2018
 *
 * USEFUL REFERENCE:
 *    -> Pthreads: https://computing.llnl.gov/tutorials/pthreads/
**********************************/

#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <fstream>
#include <pthread.h>
#include <cmath>
using namespace std;

struct thread_data{
    int thread_id;
    int thread_total;
    vector<int> *buckets;
    unsigned int numComparisons;
};

//Algorithms that are used
void insertionSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void mergeSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void quickSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void quickSortOpti(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void bucketSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void bucketSort_pthreads(vector<int> &A, unsigned int &numComparisons, clock_t &tempo, unsigned int &num_threads);

//Some auxiliary functions
void populateVectorRandom(vector<int> &A, double m, int flag);
void printVector(vector<int> A);
void printTime(clock_t &tempo);
void printComparison(unsigned int &res);
double randomNum(double m);
void exchange(int &a, int &b);
float transferTime(clock_t &tempo);
void cauculateConfidence(vector<double> &A, vector<float> &B, unsigned int &compare, unsigned int &time, double &meanTime, float &meanCompare);
void printResult(unsigned int compare, unsigned int time, float meanTime, double meanComparison);

//Sorting an Array of n integers
int main()
{
    unsigned int compare, time;
    clock_t tempo;
    unsigned int numComparisons = 0;
    double meanTime;
    float meanCompare;
    ofstream data;
    unsigned int size, flag, choice, numRun, num_threads;

//    data.open("output.txt");
//    if(data.is_open()) {
    while (true)
    {
        cout << "Please select an algorithm to run:" << endl;
        cout << "   |-->(1) insertion sort\n"
             << "   |-->(2) merge sort\n"
             << "   |-->(3) quick sort\n"
             << "   |-->(4) optimized quick sort\n"
             << "   |-->(5) bucket sort\n"
             << "   |-->(6) parallel bucket sort(implemented in pthreads)\n" << endl;
        cin >> choice;
        if(choice == 6) {
            cout << "Please input the threads number, it can be <int> or <int[]>" << endl;
            cin >> num_threads;
        }
        cout << "Please input the size of the array:" << endl;
        cin >> size;
        cout << "Please input the number: Which kind of initial array do you want?" << endl;
        cout << "   |-->(1) random array\n"
             << "   |-->(2) sorted array\n"
             << "   |-->(3) reverse array\n"
             << "   |-->(4) numerous duplications\n" << endl;
        cin >> flag;
        cout << "How many times, N, do you want to run?" << endl;
        cin >> numRun;
        cout << "Initializing....Running...." << endl;

        vector<int> A(size, 0);
        vector<double> numComparisons_data(numRun, 0);
        vector<float> time_data(numRun, 0);

        //insertion sort
        if (choice == 1)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);
                //cout << "The unsorted random Array:" << endl;
                //printVector(A1);

                //Using InsertionSort
                insertionSort(A, numComparisons, tempo);
                //                cout << "The sorted array using insertion sort:" << endl;
                //                printVector(A1);
                //                cout << "The insertion sort time cost: " << endl;
                //                printTime(tempo);
                //                cout << "The insertion sort comparison times: " << endl;
                //                printComparison(numComparisons);
                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }
        //Using MergeSort
        if (choice == 2)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);
                mergeSort(A, numComparisons, tempo);
                //                    cout << "The sorted array using merge sort:" << endl;
                //                    printVector(A2);
                //                    cout << "The merge sort time cost: " << endl;
                //                    printTime(tempo);
                //                    cout << "The merge sort comparison times: " << endl;
                //                    printComparison(numComparisons);
                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }

        //Using QuickSort
        if (choice == 3)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);
                quickSort(A, numComparisons, tempo);
                //                    cout << "The sorted array using basic quick sort:" << endl;
                //                    printVector(A3);
                //                    cout << "The quick sort time cost: " << endl;
                //                    printTime(tempo);
                //                    cout << "The quick sort comparison times: " << endl;
                //                    printComparison(numComparisons);
                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }

        //Using optimized QuickSort
        if (choice == 4)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);
                quickSortOpti(A, numComparisons, tempo);
                //                    cout << "The sorted array using optimized quick sort:" << endl;
                //                    printVector(A4);
                //                    cout << "The optimized quick sort time cost: " << endl;
                //                    printTime(tempo);
                //                    cout << "The optimized quick sort comparison times: " << endl;
                //                    printComparison(numComparisons);
                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }

        //Using BucketSort
        if (choice == 5)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);
                cout << "********************************************\n"
                     << endl;
                cout << "The unsorted array using BucketSort:\n" << endl;
                printVector(A);

                bucketSort(A, numComparisons, tempo);

                cout << "********************************************\n"
                     << endl;
                cout << "The sorted array using BucketSort:\n" << endl;
                printVector(A);

                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }

        //Using pthreads parallel BucketSort
        if (choice == 6)
        {
            for (int i = 0; i < numRun; i++)
            {
                populateVectorRandom(A, size, flag);

                populateVectorRandom(A, size, flag);
                cout << "********************************************\n"
                     << endl;
                cout << "The unsorted array using BucketSort:" << endl;
                printVector(A);

                bucketSort_pthreads(A, numComparisons, tempo, num_threads);

                cout << "********************************************\n"
                     << endl;
                cout << "The sorted array using BucketSort:" << endl;
                printVector(A);

                numComparisons_data[i] = numComparisons;
                time_data[i] = transferTime(tempo);
                numComparisons = 0;
            }
            cauculateConfidence(numComparisons_data, time_data, compare, time, meanTime, meanCompare);
            printResult(compare, time, meanTime, meanCompare);
        }
    }
    return 0;
}

//*****************************************************
//This section is an implement for InsertionSort.
//*****************************************************

void insertionSortImple(vector<int> &A, double start, double end, unsigned int &numComparisons)
{
    for (int i = start + 1; i < end; i++)
    {
        int k = A[i], j;
        for (j = i - 1; j >= 0 && A[j] > k; j--, ++numComparisons)
        {
            A[j + 1] = A[j];
        }
        ++numComparisons;
        A[j + 1] = k;
    }
}

void insertionSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo)
{
    double size = A.size();
    clock_t start, end;
    start = clock();
    insertionSortImple(A, 0, size, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for MergeSort
//*****************************************************

void mergeArray(vector<int> &A, int start, int mid, int end, int temp[], unsigned int &numComparisons)
{
    int i = start;
    int n = mid;
    int j = mid + 1;
    int m = end;
    int k = 0;
    while (i <= n && j <= m)
    {
        if (A[i] <= A[j])
        {
            ++numComparisons;
            temp[k++] = A[i++];
        }
        else
        {
            ++numComparisons;
            temp[k++] = A[j++];
        }
    }

    while (i <= n)
    {
        temp[k++] = A[i++];
    }

    while (j <= m)
    {
        temp[k++] = A[j++];
    }

    for (i = 0; i < k; i++)
    {
        A[start + i] = temp[i];
    }
}

void mergeSortSection(vector<int> &A, int start, int end, int temp[], unsigned int &numComparisons)
{
    if (start < end)
    {
        int mid = (start + end) / 2;
        mergeSortSection(A, start, mid, temp, numComparisons);
        mergeSortSection(A, mid + 1, end, temp, numComparisons);
        mergeArray(A, start, mid, end, temp, numComparisons);
    }
}

void mergeSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo)
{
    unsigned int size = A.size();
    int *temp = new int[size];
    clock_t start, end;
    start = clock();
    mergeSortSection(A, 0, size - 1, temp, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for normal quicksort
//*****************************************************

int partition(vector<int> &A, int start, int end, unsigned int &numComparisons)
{
    int key = A[end];
    int i = start - 1;

    for (int j = start; j <= end - 1; j++)
    {
        if (A[j] <= key)
        {
            ++numComparisons;
            i++;
            exchange(A[i], A[j]);
        }
    }
    exchange(A[i + 1], A[end]);
    return (i + 1);
}

void quickSortImplement(vector<int> &A, int start, int end, unsigned int &numComparisons)
{
    if (start < end)
    {
        int mid = partition(A, start, end, numComparisons);
        quickSortImplement(A, start, mid - 1, numComparisons);
        quickSortImplement(A, mid + 1, end, numComparisons);
    }
}

void quickSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo)
{
    double size = A.size();
    clock_t start, end;
    start = clock();
    quickSortImplement(A, 0, size - 1, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for optimized quicksort
//*****************************************************
int SelectPivotMedianOfThree(vector<int> &A, int low, int high, unsigned int &numComparisons)
{
    int mid = low + ((high - low) >> 1);

    //median-of-three
    if (A[mid] > A[high])
    {
        ++numComparisons;
        exchange(A[mid], A[high]);
    }
    if (A[low] > A[high])
    {
        ++numComparisons;
        exchange(A[low], A[high]);
    }
    if (A[mid] > A[low])
    {
        ++numComparisons;
        exchange(A[mid], A[low]);
    }

    return A[low];
}

int partitionOpti(vector<int> &A, int start, int end, unsigned int &numComparisons)
{
    int key = SelectPivotMedianOfThree(A, start, end, numComparisons);
    int i = start + 1, j = end;
    while (i <= j)
    {
        while (i <= end && A[i] < key)
        {
            ++numComparisons;
            i++;
        }
        while (j >= start + 1 && A[j] > key)
        {
            ++numComparisons;
            j--;
        }
        if (i > j)
            break;
        exchange(A[i], A[j]);
        i++;
        j--;
    }
    exchange(A[start], A[j]);
    return j;
}

//This section is to use tail recursion, O(n)->O(logn)
void quickSortImplementOpti(vector<int> &A, int start, int end, unsigned int &numComparisons)
{
    while (start < end)
    {
        if (end - start <= 10)
        {
            insertionSortImple(A, start, end, numComparisons);
        }
        int m = partitionOpti(A, start, end, numComparisons);
        quickSortImplementOpti(A, start, m - 1, numComparisons);
        //quickSortImplement(A, m + 1, end, numComparisons);
        start = m + 1;
    }
}

//This section is main body of enhanced quick sort
void quickSortOpti(vector<int> &A, unsigned int &numComparisons, clock_t &tempo)
{
    double size = A.size();
    //Random the pivot and convert it with the first element in the Aay
    //    double res = randomNum(size);
    //    int temp;
    //    temp = A[res];
    //    A[res] = A[0];
    //    A[0] = temp;

    clock_t start, end;
    start = clock();
    quickSortImplementOpti(A, 0, size - 1, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for BucketSort.
//*****************************************************

/**
 * Implementation of bucket sort
 */
void bucketSortImple(std::vector<int> &A, unsigned int &numComparisons)
{
    unsigned int numOfElements = A.size();

    // define num of bucket, this can choose different size
    unsigned int bucketNum = std::ceil(numOfElements / 2.0);
    int minValue = A[0], maxValue = A[0];

    for (int i = 0; i < numOfElements; i++)
    {
        if (A[i] > maxValue)
        {
            maxValue = A[i];
        }
        if (A[i] < minValue)
        {
            minValue = A[i];
        }
    }

    // range of each bucket that will persist for the element in the array
    double rangePerBucket = (maxValue - minValue + 1) / (double)bucketNum;

    std::vector<int> *buckets = new std::vector<int>[bucketNum];

    for (int i = 0; i < bucketNum; i++)
    {
        buckets[i] = std::vector<int>();
    }

    // push each element into bucket
    for (int i = 0; i < numOfElements; i++)
    {

        int bucketIndex = std::floor((A[i] - minValue) / rangePerBucket);

        buckets[bucketIndex].push_back(A[i]);
    }

    // sort each bucket's element in specific sorting algorithm, here uses Insertion sort
    for (int i = 0; i < bucketNum; i++)
    {
        if (buckets[i].size() > 0)
        {
//            cout << "Each elements in bucket[" << i << "]: " << endl;
            insertionSortImple(buckets[i], 0, buckets[i].size(), numComparisons);
//            printVector(buckets[i]);
        }
    }

    // reduce the elements back to original array
    int index = 0;
    for (int i = 0; i < bucketNum; i++)
    {
        for (int j = 0; j < buckets[i].size(); j++)
        {
            A[index++] = buckets[i][j];
        }
    }
}

void bucketSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo)
{
    unsigned int size = A.size();
    clock_t start, end;
    start = clock();
    bucketSortImple(A, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for parallel BucketSort using pthreads.
//*****************************************************

/**
 * insertion sort wrapper with pthreads
 */
void* thread_insertion_sort(void* threadarg) {

    struct thread_data *local_data = (struct thread_data *) threadarg;
    int partition = local_data->thread_id;
    int thread_total = local_data->thread_total;
    std::vector<int>* buckets = local_data->buckets;

    // calculate the partition range
    int range = buckets->size() / thread_total;
    if(range * thread_total < buckets->size()) {
        range += 1;
    }

    if(partition == thread_total - 1) {
        for(int i = partition * range; i < buckets->size(); i++) {
            insertionSortImple(*buckets, 0, buckets[i].size(), local_data->numComparisons);
        }
    } else {
        for(int i = partition * range; i < (partition + 1) * range; i++) {
            insertionSortImple(*buckets, 0, buckets[i].size(), local_data->numComparisons);
        }
    }

    pthread_exit(NULL);
}

/**
 * Worker threads creator
 */

void threads_factory(std::vector<int> &buckets, int num_threads, unsigned int &numComparisons) {
    pthread_t threads[num_threads];
    struct thread_data thread_data_array[num_threads];

    // create threads and start work
    cout << "creating "<< num_threads << " threads........." << endl;
    for(int tid = 0; tid < num_threads; tid++) {
        thread_data_array[tid].thread_id = tid;
        thread_data_array[tid].buckets = &buckets;
        thread_data_array[tid].thread_total = num_threads;
        thread_data_array[tid].numComparisons = 0;

        int rc = pthread_create(&threads[tid], NULL, thread_insertion_sort, (void *) &thread_data_array[tid]);

        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    // wait for all threads to complete their work
    for(int tid = 0; tid < num_threads; tid++) {
        int rc = pthread_join(threads[tid], NULL);
        if (rc)
        {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            return;
        }
    }

    for(int tid = 0; tid < num_threads; tid++) {
        numComparisons += thread_data_array[tid].numComparisons;
    }
}

/**
 * Implementation of parallel bucket sort using pthreads
 */
void bucketSort_pthreads_imple(std::vector<int> &A, unsigned int &numComparisons, unsigned int &num_threads)
{
    unsigned int numOfElements = A.size();

    // define num of bucket, this can choose different size
    unsigned int bucketNum = std::ceil(numOfElements / 2.0);
    int minValue = A[0], maxValue = A[0];

    for (int i = 0; i < numOfElements; i++)
    {
        if (A[i] > maxValue)
        {
            maxValue = A[i];
        }
        if (A[i] < minValue)
        {
            minValue = A[i];
        }
    }

    // range of each bucket that will persist for the element in the array
    double rangePerBucket = (maxValue - minValue + 1) / (double)bucketNum;

    std::vector<int> *buckets = new std::vector<int>[bucketNum];

    for (int i = 0; i < bucketNum; i++)
    {
        buckets[i] = std::vector<int>();
    }

    // push each element into bucket
    for (int i = 0; i < numOfElements; i++)
    {

        int bucketIndex = static_cast<int>(std::floor((A[i] - minValue) / rangePerBucket));

        buckets[bucketIndex].push_back(A[i]);
    }

    // start thread workers
    threads_factory(*buckets, num_threads, numComparisons);

    // reduce the elements back to original array
    int index = 0;
    for (int i = 0; i < bucketNum; i++)
    {
        for (int j = 0; j < buckets[i].size(); j++)
        {
            A[index++] = buckets[i][j];
        }
    }
}

void bucketSort_pthreads(vector<int> &A, unsigned int &numComparisons, clock_t &tempo, unsigned int &num_threads)
{
    clock_t start, end;
    start = clock();
    bucketSort_pthreads_imple(A, numComparisons, num_threads);
    end = clock();
    tempo = end - start;
}

/*****************************************************
 * Helper functions
 *****************************************************/

//This section is to populate the Array with
//integer numbers selected randomly and uniformly in the range [0, m].
void populateVectorRandom(vector<int> &A, double m, int flag)
{
    //Generate random sequence
    random_device rd;
    default_random_engine e(rd());
    uniform_int_distribution<> range1(0, 10 * (m + 1));
    uniform_int_distribution<> range2(0, (int)A.size() / 50);
    uniform_int_distribution<> range3(0, 10);

    //Pupulate random numbers into vectors
    if (flag == 1)
    {
        for (double i = 0; i < A.size(); i++)
        {
            A[i] = range1(e);
        }
    }
    //pupulate sorted array
    else if (flag == 2)
    {
        for (double i = 0; i < A.size(); i++)
        {
            A[i] = i + range3(e) + 1;
        }
    }
    //populate reversed array
    else if (flag == 3)
    {
        for (double i = 0; i < A.size(); i++)
        {
            A[i] = A.size() * 2 - range3(e) - i;
        }
    }
    //populate few unique array
    else if (flag == 4)
    {
        for (double i = 0; i < A.size(); i++)
        {
            A[i] = range3(e);
        }
    }
}

//This section is to print an Aay out.
void printVector(vector<int> A)
{
    int n = 0;
    for (int &i : A)
    {
        // Set a threshold to print
        if (n >= 20) {
            cout << "only showing "<< n << " elements, omit remaining......";
            break;
        }
        cout << i << "\t";

        //Each line with 5 elements.
        if (++n % 5 == 0)
            cout << "\n";
    }
    cout << endl;
}

//Print the time
void printTime(clock_t &tempo)
{
    cout << (float)(tempo) / CLOCKS_PER_SEC << endl;
}

//Print the comparison
void printComparison(unsigned int &res)
{
    cout << res << endl;
}

//Transfer time
float transferTime(clock_t &tempo)
{
    return (float)tempo / CLOCKS_PER_SEC;
}

//generate random number
double randomNum(double m)
{
    //Generate random num
    random_device rd;
    default_random_engine e(rd());
    uniform_int_distribution<> range(0, m);
    return range(e);
}

//exchange a with b
void exchange(int &a, int &b)
{
    int temp;
    temp = a;
    a = b;
    b = temp;
}

void cauculateConfidence(vector<double> &A, vector<float> &B, unsigned int &compare, unsigned int &time, double &meanTime, float &meanCompare)
{
    double x1 = 0, s1 = 0;
    float x2 = 0.0, s2 = 0.0;
    x1 = accumulate(begin(A), end(A), 0.0);
    x1 = x1 / A.size();
    x2 = static_cast<float>(accumulate(begin(B), end(B), 0.0));
    x2 = x2 / B.size();
    for (int i = 0; i < 50; i++)
    {
        s1 += pow(A[i] - x1, 2);
        s2 += pow(B[i] - x2, 2);
    }
    s1 = s1 / (A.size() - 1);
    s2 = s2 / (B.size() - 1);

    meanCompare = static_cast<float>(x1);
    meanTime = x2;

    compare = (unsigned int)(1536.64 * s1 / pow(x1, 2));
    time = (unsigned int)(1536.64 * s2 / pow(x2, 2));
    //cout << time << endl;

    //    for(int i = 0; i < 50; i++) {
    //        x1 += A[i];
    //        x2 += B[i];
    //    }
    //
    //    x2 = x2 / 50;
}

void printResult(unsigned int compare, unsigned int time, float meanTime, double meanComparison)
{
//    cout << "The n for Comparisons: " << endl;
//    cout << compare << endl;
//    cout << "The n for Time: " << endl;
//    cout << time << endl;
    cout << "The mean value of time : " << endl;
    cout << meanTime << endl;
    cout << "The mean value of comparisons: " << endl;
    cout << meanComparison << endl;
}

//                    Print instruction of this program.
//                    cout << "Attention that n is an integer number, n > 1 and m is a real number, m > 0\n" << endl;
//                    cout << "Please input the size of the Aay:" << endl;
//                    cin >> n;
//                    while (n <= 1)
//                    {
//                        cout << "Wrong input, the size of Aay must be more than 1" << endl;
//                        cout << "Please input the size of the Aay:" << endl;
//                        cin >> n;
//                    }
//
//                    cout << "Please input the range:" << endl;
//                    cin >> m;
//                    while (m <= 1)
//                    {
//                        cout << "Wrong input, the range must be more than 1" << endl;
//                        cout << "Please input the range:" << endl;
//                        cin >> m;
//                    }
