// Include your C header files here

#include "stdio.h"
#include "stdlib.h"
#include "pthread.h"
#include "pth_msort.h"

struct ArrayStruct
{
    int *arr;
    int left;
    int right;
};

struct MergeStruct
{
    int thread_ID;
    int *arr1;
    int *arr2;
    int *merged_arr;
    int length;
};

void merge_sort(int arr[], int l, int r)
{
    if (l < r)
    {
        int m = l + (r - l) / 2;

        merge_sort(arr, l, m);
        merge_sort(arr, m + 1, r);

        // merge
        int i, j, k;
        int n1 = m - l + 1;
        int n2 = r - m;

        int *tempL = (int *)malloc(sizeof(int) * n1);
        int *tempR = (int *)malloc(sizeof(int) * n2);

        for (i = 0; i < n1; i++)
            tempL[i] = arr[l + i];
        for (j = 0; j < n2; j++)
            tempR[j] = arr[m + 1 + j];

        i = 0;
        j = 0;
        k = l;
        while (i < n1 && j < n2)
            if (tempL[i] <= tempR[j])
                arr[k++] = tempL[i++];
            else
                arr[k++] = tempR[j++];

        while (i < n1)
            arr[k++] = tempL[i++];

        while (j < n2)
            arr[k++] = tempR[j++];

        free(tempL);
        free(tempR);
    }
}

void *MergeSort(void *Array)
{
    int *arr = ((struct ArrayStruct *)Array)->arr;
    int left = ((struct ArrayStruct *)Array)->left;
    int right = ((struct ArrayStruct *)Array)->right;
    merge_sort(arr, left, right);
}

void *SerialMerge(void *StructedArray)
{
    int *arr1 = ((struct MergeStruct *)StructedArray)->arr1;
    int *arr2 = ((struct MergeStruct *)StructedArray)->arr2;
    int *merged_arr = ((struct MergeStruct *)StructedArray)->merged_arr;
    int length = ((struct MergeStruct *)StructedArray)->length;

    int i = 0, j = 0, k = 0;
    while (i < length && j < length)
        if (arr1[i] < arr2[j])
            merged_arr[k++] = arr1[i++];
        else
            merged_arr[k++] = arr2[j++];
    while (i < length)
        merged_arr[k++] = arr1[i++];
    while (j < length)
        merged_arr[k++] = arr2[j++];
}

int BinarySearch1(int arr[], int l, int r, int x)
{
    int mid;
    while (r > l)
    {
        mid = l + ((r - l) / 2);
        if (arr[mid] < x)
            l = mid + 1;
        else
            r = mid;
    }
    if (x > arr[r])
        return r + 1;
    return r;
}

int BinarySearch2(int arr[], int l, int r, int x)
{
    int mid;
    while (r > l)
    {
        mid = l + ((r - l) / 2);
        if (arr[mid] > x)
            r = mid;
        else
            l = mid + 1;
    }
    if (x >= arr[r])
        return r + 1;
    return r;
}

void *ParallelMerge(void *StructedArray)
{
    int thread_ID = ((struct MergeStruct *)StructedArray)->thread_ID;
    int *Arr1 = ((struct MergeStruct *)StructedArray)->arr1;
    int *Arr2 = ((struct MergeStruct *)StructedArray)->arr2;
    int *merged_arr = ((struct MergeStruct *)StructedArray)->merged_arr;
    int length = ((struct MergeStruct *)StructedArray)->length;

    int arr1_left = thread_ID * length / 8;
    int arr1_right = (thread_ID + 1) * (length / 8) - 1;

    int i, j = 0, k = 0;

    for (i = arr1_left; i <= arr1_right; i++)
    {
        j = BinarySearch1(Arr2, j, (length / 2 - 1), Arr1[i]);
        merged_arr[j + i] = Arr1[i];

        k = BinarySearch2(Arr1, k, (length / 2 - 1), Arr2[i]);
        merged_arr[k + i] = Arr2[i];
    }
}

void mergeArrays(int arr1[], int arr2[], int n1, int n2, int arr3[])
{
    int i = 0, j = 0, k = 0;

    while (i < n1 && j < n2)
        if (arr1[i] < arr2[j])
            arr3[k++] = arr1[i++];
        else
            arr3[k++] = arr2[j++];

    while (i < n1)
        arr3[k++] = arr1[i++];

    while (j < n2)
        arr3[k++] = arr2[j++];
}

void mergeSortParallel(const int *values, unsigned int N, int *sorted)
{
    int *arr_values = (int *)values;

    if (N >> 30)
    {
        mergeSortParallel(values, N / 2, arr_values);
        mergeSortParallel(&values[N / 2], N / 2, &arr_values[N / 2]);

        mergeArrays(&arr_values[0], &arr_values[N / 2], N / 2, N / 2, sorted);
        return;
    }

    int i = 0;
    int T = 4;

    // parallel quicksort //
    struct ArrayStruct Arr[4] = {{arr_values, 0, N / 4 - 1},
                                 {arr_values, N / 4, N / 2 - 1},
                                 {arr_values, N / 2, 3 * N / 4 - 1},
                                 {arr_values, 3 * N / 4, N - 1}};

    pthread_t *handles = (pthread_t *)malloc(T * sizeof(pthread_t));

    for (i = 0; i < T; i++)
        pthread_create(&handles[i], NULL, MergeSort, (void *)&Arr[i]);

    for (i = 0; i < T; i++)
        pthread_join(handles[i], NULL);

    // 2 middle merge modules //
    int *MergeTemp = (int *)malloc(sizeof(int) * N);
    struct MergeStruct merged[2] = {{0, &arr_values[0], &arr_values[N / 4], &MergeTemp[0], N / 4},
                                    {1, &arr_values[N / 2], &arr_values[3 * N / 4], &MergeTemp[N / 2], N / 4}};

    pthread_create(&handles[0], NULL, SerialMerge, (void *)&merged[0]);
    pthread_create(&handles[1], NULL, SerialMerge, (void *)&merged[1]);
    pthread_join(handles[0], NULL);
    pthread_join(handles[1], NULL);

    // parallel merge sort //
    struct MergeStruct SortedMerge[4] = {{0, MergeTemp, &MergeTemp[N / 2], sorted, N},
                                         {1, MergeTemp, &MergeTemp[N / 2], sorted, N},
                                         {2, MergeTemp, &MergeTemp[N / 2], sorted, N},
                                         {3, MergeTemp, &MergeTemp[N / 2], sorted, N}};

    for (i = 0; i < T; i++)
        pthread_create(&handles[i], NULL, ParallelMerge, (void *)&SortedMerge[i]);

    for (i = 0; i < T; i++)
        pthread_join(handles[i], NULL);

    free(handles);
    free(MergeTemp);
}