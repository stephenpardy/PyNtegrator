#include <stdio.h>
#include <stdlib.h>

int test(int *arr){
    int i;
    if (arr){
        for (i=0; i<10; i++) printf("%d\n", arr[i]);
    }
}


int main(void){
    int arr[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    test(arr);
    test(NULL);
    test(&arr[0]);
}