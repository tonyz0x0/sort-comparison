//
//  main.cpp
//  AlgorithmComparison
//
//  Created by Kejie Zhang on 2017/11/22.
//  Copyright © 2017 Kejie Zhang. All rights reserved.
//

#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <fstream>
using namespace std;

//All the functions that are used

//Four algorithms that are used
void insertionSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void mergeSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void quickSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);
void quickSortOpti(vector<int> &A, unsigned int &numComparisons, clock_t &tempo);


//Some auxiliary functions
void populateVectorRandom(vector<int> &A, double m, int flag);
void printVector (vector<int> A);
void printTime (clock_t &tempo);
void printComparison(unsigned int &res);
double randomNum(double m);
void exchange(int &a, int &b);
float transferTime(clock_t &tempo);
void cauculateConfidence(vector<double> &A, vector<float> &B, unsigned int &compare, unsigned int &time, double &meanTime, float &meanCompare);
void printResult(unsigned int compare, unsigned int time, float meanTime, double meanComparison);

//Sorting an Aay of n integers
int main()
{
    vector<double> A7(50, 0);
    vector<float> A8(50, 0);
    unsigned int compare, time;
    clock_t tempo1, tempo2, tempo3, tempo4;
    unsigned int numComparisons1 = 0, numComparisons2 = 0, numComparisons3 = 0, numComparisons4 = 0;
    double meanTime;
    float meanCompare;
    ofstream data;
    unsigned int size, flag, choice, numRun;
    
    
    data.open("output.txt");
    //if(data.is_open()) {
        while(true) {
            cout << "Please select an algorithm to run:" << endl;
            cout << "(1)insertion sort; (2)merge sort; (3)quick sort; (4)optimized quick sort" << endl;
            cin >> choice;
            cout << "Please input the size of the array:" << endl;
            cin >> size;
            cout << "Please input the number: Which kind of initial array do you want?" << endl;
            cout << "(1)random array; (2)sorted array; (3)reverse array; (4)numerous duplications" << endl;
            cin >> flag;
            cout << "How many times, N, do you want to run?" << endl;
            cin >> numRun;

        vector<int> A(size, 0);
            
            //insertion sort
            if(choice == 1){
                for(int i = 0; i < numRun; i++) {
                        populateVectorRandom(A, size, flag);
                        //cout << "The unsorted random Array:" << endl;
                        //printVector(A1);

                        //Using InsertionSort
                        insertionSort(A, numComparisons1, tempo1);
        //                cout << "The sorted array using insertion sort:" << endl;
        //                printVector(A1);
        //                cout << "The insertion sort time cost: " << endl;
        //                printTime(tempo1);
        //                cout << "The insertion sort comparison times: " << endl;
        //                printComparison(numComparisons1);
                        A7[i] = numComparisons1;
                        A8[i] = transferTime(tempo1);
                    numComparisons1 = 0;
                      }
                cauculateConfidence(A7, A8, compare, time, meanTime, meanCompare);
                printResult(compare, time, meanTime, meanCompare);
            }
                //Using MergeSort
            if(choice == 2){
                for(int i = 0; i < numRun; i++) {
                    populateVectorRandom(A, size, flag);
                    mergeSort(A,numComparisons2, tempo2);
//                    cout << "The sorted array using merge sort:" << endl;
//                    printVector(A2);
//                    cout << "The merge sort time cost: " << endl;
//                    printTime(tempo2);
//                    cout << "The merge sort comparison times: " << endl;
//                    printComparison(numComparisons2);
                    A7[i] = numComparisons2;
                    A8[i] = transferTime(tempo2);
                    numComparisons2 = 0;
                }
                cauculateConfidence(A7, A8, compare, time, meanTime, meanCompare);
                printResult(compare, time, meanTime, meanCompare);
            }
            
            
                //Using QuickSort
            if(choice == 3){
                for(int i = 0; i < numRun; i++) {
                    populateVectorRandom(A, size, flag);
                    quickSort(A, numComparisons3, tempo3);
//                    cout << "The sorted array using basic quick sort:" << endl;
//                    printVector(A3);
//                    cout << "The quick sort time cost: " << endl;
//                    printTime(tempo3);
//                    cout << "The quick sort comparison times: " << endl;
//                    printComparison(numComparisons3);
                    A7[i] = numComparisons3;
                    A8[i] = transferTime(tempo3);
                    numComparisons3 = 0;
                }
                cauculateConfidence(A7, A8, compare, time, meanTime, meanCompare);
                printResult(compare, time, meanTime, meanCompare);
                
            }
                
                //Using optimized QuickSort
            if(choice == 4){
                for(int i = 0; i < numRun; i++) {
                    populateVectorRandom(A, size, flag);
                    quickSortOpti(A, numComparisons4, tempo4);
//                    cout << "The sorted array using optimized quick sort:" << endl;
//                    printVector(A4);
//                    cout << "The optimized quick sort time cost: " << endl;
//                    printTime(tempo4);
//                    cout << "The optimized quick sort comparison times: " << endl;
//                    printComparison(numComparisons4);
                    A7[i] = numComparisons4;
                    A8[i] = transferTime(tempo4);
                    numComparisons4 = 0;
                }
                cauculateConfidence(A7, A8, compare, time, meanTime, meanCompare);
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
    for(int i = start + 1; i < end; i++)
    {
        int k = A[i], j;
        for(j = i - 1; j >= 0 && A[j] > k; j--, ++numComparisons)
        {
            A[j + 1] = A[j];
        }
        ++numComparisons;
        A[j + 1] = k;
    }
}

void insertionSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo){
    double size = A.size();
    clock_t start, end;
    start = clock();
    insertionSortImple(A,0, size, numComparisons);
    end = clock();
    tempo = end - start;
}

//*****************************************************
//This section is an implement for MergeSort
//*****************************************************

void mergeArray(vector<int> &A, int start, int mid, int end, int temp[], unsigned int &numComparisons) {
    int i = start;
    int n = mid;
    int j = mid + 1;
    int m = end;
    int k = 0;
    while(i <= n && j <= m) {
        if(A[i] <= A[j]) {
            ++numComparisons;
            temp[k++] = A[i++];
        } else {
            ++numComparisons;
            temp[k++] = A[j++];
        }
    }
    
    while(i <= n) {
        temp[k++] = A[i++];
    }
    
    while(j <= m) {
        temp[k++] = A[j++];
    }
    
    for(i = 0; i < k; i++) {
        A[start + i] = temp[i];
    }
}

void mergeSortSection(vector<int> &A, int start, int end, int temp[], unsigned int &numComparisons) {
    if(start < end) {
        int mid = (start + end) / 2;
        mergeSortSection(A, start, mid, temp, numComparisons);
        mergeSortSection(A, mid + 1, end, temp, numComparisons);
        mergeArray(A, start, mid, end, temp, numComparisons);
    }
}

void mergeSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo) {
    double size = A.size();
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

int partition (vector<int> &A, int start, int end, unsigned int &numComparisons) {
    int key = A[end];
    int i = start - 1;
    
    for(int j = start; j <= end - 1; j ++) {
        if(A[j] <= key) {
            ++numComparisons;
            i++;
            exchange(A[i], A[j]);
        }
    }
    exchange(A[i + 1], A[end]);
    return (i + 1);
}

void quickSortImplement(vector<int> &A, int start, int end, unsigned int &numComparisons) {
    if(start < end) {
        int mid = partition(A, start, end, numComparisons);
        quickSortImplement(A, start, mid - 1, numComparisons);
        quickSortImplement(A, mid + 1, end, numComparisons);
    }
}

void quickSort(vector<int> &A, unsigned int &numComparisons, clock_t &tempo){
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
int SelectPivotMedianOfThree(vector<int> &A,int low,int high, unsigned int &numComparisons)
{
    int mid = low + ((high - low) >> 1);
    
    //median-of-three
    if (A[mid] > A[high])
    {
        ++numComparisons;
        exchange(A[mid],A[high]);
    }
    if (A[low] > A[high])
    {
        ++numComparisons;
        exchange(A[low],A[high]);
    }
    if (A[mid] > A[low])
    {
        ++numComparisons;
        exchange(A[mid],A[low]);
    }
    
    return A[low];
}


int partitionOpti(vector<int> &A, int start, int end, unsigned int &numComparisons){
    int key =SelectPivotMedianOfThree(A, start, end, numComparisons);
    int i = start + 1, j = end;
    while (i <= j) {
        while(i <= end && A[i] < key) {
            ++numComparisons;
            i++;
        }
        while(j >= start + 1 && A[j] > key) {
            ++numComparisons;
            j--;
        }
        if(i > j) break;
        exchange(A[i], A[j]);
        i++;
        j--;
    }
    exchange(A[start], A[j]);
    return j;
}


//This section is to use tail recursion, O(n)->O(logn)
void quickSortImplementOpti(vector<int> &A, int start, int end, unsigned int &numComparisons) {
    while(start < end)
    {
        if(end - start <= 10) {
            insertionSortImple(A, start, end, numComparisons);
        }
        int m = partitionOpti(A, start, end, numComparisons);
        quickSortImplementOpti(A, start, m - 1, numComparisons);
        //quickSortImplement(A, m + 1, end, numComparisons);
        start = m + 1;
    }
}

//This section is main body of enhanced quick sort
void quickSortOpti(vector<int> &A, unsigned int &numComparisons, clock_t &tempo) {
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


//This section is to populate the Array with
//integer numbers selected randomly and uniformly in the range [0, m].
void populateVectorRandom(vector<int> &A, double m, int flag)
{
    //Generate random sequence
    random_device rd;
    default_random_engine e(rd());
    uniform_int_distribution<> range1(0, 10*(m + 1));
    uniform_int_distribution<> range2(0, (int)A.size() / 50);
    uniform_int_distribution<> range3(0, 10);
    
    //Pupulate random numbers into vectors
    if(flag == 1) {
        for(double i = 0; i < A.size(); i++)
        {
            A[i] = range1(e);
        }
    }
    //pupulate sorted array
    else if (flag == 2) {
        for(double i = 0; i < A.size(); i++)
        {
            A[i] = i + range3(e) + 1;
        }
    }
    //populate reversed array
    else if (flag == 3) {
        for(double i = 0; i < A.size(); i++)
        {
            A[i] = A.size() * 2 - range3(e)- i;
        }
    }
    //populate few unique array
    else if (flag == 4) {
        for(double i = 0; i < A.size(); i++)
        {
            A[i] = range3(e);
        }
    }
}

//This section is to print an Aay out.
void printVector (vector<int> A)
{
    int n = 0;
    for(int &i : A)
    {
        cout << i << "\t";
        
        //Each line with 5 elements.
        if(++n % 5 == 0)
            cout << "\n";
    }
    cout << endl;
}

//Print the time
void printTime (clock_t &tempo){
    cout << (float)(tempo) / CLOCKS_PER_SEC << endl;
}

//Print the comparison
void printComparison(unsigned int &res){
    cout << res << endl;
}

//Transfer time
float transferTime(clock_t &tempo) {
    return (float)tempo/CLOCKS_PER_SEC;
}

//generate random number
double randomNum(double m) {
    //Generate random num
    random_device rd;
    default_random_engine e(rd());
    uniform_int_distribution<> range(0, m);
    return range(e);
}

//exchange a with b
void exchange(int &a, int &b){
    int temp;
    temp = a;
    a = b;
    b = temp;
}

void cauculateConfidence(vector<double> &A, vector<float> &B, unsigned int &compare, unsigned int &time, double &meanTime, float &meanCompare) {
    double x1 = 0, s1 = 0;
    float x2 = 0.0, s2 = 0.0;
    x1 = accumulate(begin(A), end(A), 0.0);
    x1 = x1 / A.size();
    x2 = accumulate(begin(B), end(B), 0.0);
    x2 = x2 / B.size();
    for(int i = 0; i < 50; i++) {
        s1 += pow(A[i] - x1, 2);
        s2 += pow(B[i] - x2, 2);
    }
    s1 = s1 / (A.size() - 1);
    s2 = s2 / (B.size() - 1);
    
    meanTime = x1;
    meanCompare = x2;
    
    compare = (int)(1536.64 * s1 / pow(x1, 2));
    time = (int)(1536.64 * s2 / pow(x2, 2));
    //cout << time << endl;
    
//    for(int i = 0; i < 50; i++) {
//        x1 += A[i];
//        x2 += B[i];
//    }
//
//    x2 = x2 / 50;
}

void printResult(unsigned int compare, unsigned int time, float meanTime, double meanComparison) {
    cout << "The n for Comparisons: " << endl;
    cout << compare << endl;
    cout << "The n for Time: " << endl;
    cout << time << endl;
    cout << "The mean value comparisons : " << endl;
    cout << meanTime << endl;
    cout << "The mean value of time: " << endl;
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
    
