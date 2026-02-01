#include <stdio.h>

void bubble_sort(int *arr, int n) {
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (arr[j] > arr[j + 1]) {
                int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

int main() {
    int data[] = {5, 2, 9, 1, 3};
    int size = sizeof(data) / sizeof(data[0]);

    bubble_sort(data, size);

    printf("冒泡排序结果：");
    for (int i = 0; i < size; ++i) {
        printf("%d ", data[i]);
    }
    printf("\n");

    return 0;
}