#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

struct Matrix{
    int startRow;
    int startColumn;
    int dimension;
    int** values;
};

int** strassen(int** m1, int** m2, int n);
int** bruteForce(int** m1, int** m2, int n);
Matrix* initMatrix(int d);
void populateMatrices(Matrix* A, Matrix* B, int d, char* inputfile);
Matrix* addMatrices(Matrix* A, Matrix* B, bool subtract);
void printMat(Matrix* matrix, int d);

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
    Matrix* A = initMatrix(d);
    Matrix* B = initMatrix(d);

    populateMatrices(A, B, d, inputfile);

    printMat(A, d);
    printMat(B, d);

    return 0;
}

Matrix* initMatrix(int d) {
    Matrix* mat = (Matrix*) malloc(sizeof(Matrix));
    mat->dimension = d;
    mat->startRow = 0;
    mat->startColumn = 0;
    int** A = (int**) malloc(sizeof(int*) * d);
    for (int i = 0; i < d; i++) {
        A[i] = (int*) malloc(d * sizeof(int));
    }
    mat->values = A;
    return mat;
}

void populateMatrices(Matrix* A, Matrix* B, int d, char* inputfile) {
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
            A->values[i][j] = el;
        }
    }
    // populate matrix B
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            fscanf(file, "%s\n", line);
            int el = strtol(line, NULL, 10);
            B->values[i][j] = el;
        }
    }
    fclose(file);
}

void printMat(Matrix* mat, int d) {
  for (int i = mat->startRow; i < d; i++) {
    for (int j = mat->startColumn; j < d; j++)
      printf("%i ", mat->values[i][j]);
    printf("\n");
  }
}

Matrix* addMatrices(Matrix* A, Matrix* B, bool subtract) {
    if (A->dimension != B->dimension) {
        printf("Error with adding matrices.\n");
        return;
    }
    int d = A->dimension;
    Matrix* res = initMatrix(d);
    int ARidx = A->startRow;
    int ACidx = A->startColumn;
    int BRidx = B->startRow;
    int BCidx = B->startColumn;
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (subtract) {
                res->values[i][j] = A->values[ARidx + i][ACidx + j] - B->values[BRidx + i][BCidx + j];
            } else {
                res->values[i][j] = A->values[ARidx + i][ACidx + j] + B->values[BRidx + i][BCidx + j];
            }
        }
    }
    return res;
}

Matrix* splitMatrices()