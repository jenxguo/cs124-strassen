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

Matrix* strassen(Matrix* m1, Matrix* m2);
Matrix* conventionalMult(Matrix* m1, Matrix* m2);
Matrix* initMatrix(int d);
Matrix* splitMatrix(int startRow, int startCol);
void populateMatrices(Matrix* A, Matrix* B, int d, char* inputfile);
void printMat(Matrix* matrix);

int N_0 = 15;

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

    printMat(A);
    printMat(B);

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

void printMat(Matrix* mat) {
    int d = mat->dimension;
    for (int i = mat->startRow; i < d; i++) {
        for (int j = mat->startColumn; j < d; j++)
        printf("%i ", mat->values[i][j]);
        printf("\n");
    }
}

Matrix* addMatrices(Matrix* A, Matrix* B, bool subtract) {
    if (A->dimension != B->dimension) {
        printf("Error with adding matrices.\n");
        return NULL;
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

// Matrix* splitMatrix(int startRow, int startCol);

// Matrix* createMatrix(int startRow, int startCol, ){

// };


Matrix* strassen(Matrix* m1, Matrix* m2){
    if (m1->dimension == N_0) {
        return conventionalMult(m1, m1);
    }

    Matrix* A;
    Matrix* B;
    Matrix* C;
    Matrix* D;
    Matrix* E;
    Matrix* F;
    Matrix* G;
    Matrix* H;

    Matrix* P1 = (Matrix*) malloc(sizeof(Matrix));
    P1 = strassen(A, addMatrices(F, H, true));

    Matrix* P2 = (Matrix*) malloc(sizeof(Matrix));
    P2 = strassen(addMatrices(A, B, false), H);

    Matrix* P3 = (Matrix*) malloc(sizeof(Matrix));
    P3 = strassen(addMatrices(C, D, false), E);

    Matrix* P4 = (Matrix*) malloc(sizeof(Matrix));
    P4 = strassen(D, addMatrices(G, E, true));

    Matrix* P5 = (Matrix*) malloc(sizeof(Matrix));
    P5 = strassen(addMatrices(A, D, false), addMatrices(E, H, false));

    Matrix* P6 = (Matrix*) malloc(sizeof(Matrix));
    P6 = strassen(addMatrices(B, D, false), addMatrices(G, H, false));

    Matrix* P7 = (Matrix*) malloc(sizeof(Matrix));
    P7 = strassen(addMatrices(C, A, true), addMatrices(E, F, false));
};


Matrix* conventionalMult(Matrix* m1, Matrix* m2){
    int d = m1->dimension;
    int ARidx = m1->startRow;
    int ACidx = m1->startColumn;
    int BRidx = m2->startRow;
    int BCidx = m2->startColumn;
    Matrix* res = initMatrix(d);
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            res->values[i][j] = 0;
            for (int k = 0; k < d; k++) {
                res->values[i][j] += m1->values[ARidx + i][ACidx + k] * m2->values[BRidx + k][BCidx + j];
            }
        }
    }
    return res;
};
