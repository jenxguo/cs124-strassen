#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

int** strassen(int** m1, int** m2, int n);
int** bruteForce(int** m1, int** m2, int n);
void populateMatrices(int** A, int** B, int d, char* inputfile);
void printMat(int** matrix, int d);

struct Matrix{
    int startRow;
    int startColumn;
    int dimension;
    int** values;
};

int main(int argc, char *argv[]) {

    //checks that there are command line args
    if (argc != 4)
    {
        cout << "Usage: ./randmst 0 dimension inputfile\n";
        return 1;
    }

    int d = (int) strtol(argv[2], NULL, 10);
    char *inputfile = argv[3];

    // https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
    // if dimension is power of 2
    int** A = (int**) malloc(sizeof(int*) * d);
    for (int i = 0; i < d; i++) {
        A[i] = (int*) malloc(d * sizeof(int));
    }
    int** B = (int**) malloc(sizeof(int*) * d);
    for (int i = 0; i < d; i++) {
        B[i] = (int*) malloc(d * sizeof(int));
    }

    populateMatrices(A, B, d, inputfile);

    printMat(A, d);
    printMat(B, d);

    return 0;
}

void populateMatrices(int** A, int** B, int d, char* inputfile) {
    FILE* file = fopen(inputfile, "r");
    if (file == NULL) {
        printf("File could not be opened.\n");
        return;
    }
    char line[10];
    // populate matrix A
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            fscanf(file, "%s\n", line);
            int el = strtol(line, NULL, 10);
            A[i][j] = el;
        }
    }
    // populate matrix B
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            fscanf(file, "%s\n", line);
            int el = strtol(line, NULL, 10);
            B[i][j] = el;
        }
    }
    fclose(file);
}

void printMat(int** matrix, int d) {
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j)
      printf("%i ", matrix[i][j]);
    printf("\n");
  }
}
