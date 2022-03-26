#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <tuple>
using namespace std;


struct Matrix{
    int startRow;
    int startColumn;
    int dimension;
    int** values;
};

int calcPadding(int d);
Matrix* strassen(Matrix* m1, Matrix* m2);
Matrix* conventionalMult(Matrix* m1, Matrix* m2);
Matrix** splitMatrices(Matrix* original);
void combine(Matrix* Product, Matrix* TL, Matrix* TR, Matrix* BL, Matrix* BR);
Matrix* addMatrices(Matrix* A, Matrix* B, bool subtract);
void populateMatrices(int** A, int** B, int d, char* inputfile);
void printMat(int** matrix, int d);
Matrix* initMatrix(int d);
void freeMatrix(Matrix* mat);

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

    int pad = calcPadding(d);

    Matrix* A = initMatrix(pad);
    Matrix* B = initMatrix(pad);

    populateMatrices(A, B, d, inputfile);

    printMat(A);
    printMat(B);

    return 0;
}

int calcPadding(int d) {
    if (d <= N_0) {
        return d;
    }
    int n = d;
    while (n > N_0) {
        n = (n / 2) + ((n % 2) != 0) ;
    }
    while (n < d) {
        n = n*2;
    }
    return n;
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


Matrix** splitMatrices(Matrix* original){
    int n = original->dimension;
    int newDim = n / 2;

    Matrix* A = (Matrix*) malloc(sizeof(Matrix));
    A->dimension = newDim;
    A->startRow = original->startRow;
    A->startColumn = original->startColumn;
    A->values = original->values;

    Matrix* B = (Matrix*) malloc(sizeof(Matrix));
    B->dimension = newDim;
    B->startRow = original->startRow;
    B->startColumn = original->startColumn + newDim;
    B->values = original->values;

    Matrix* C = (Matrix*) malloc(sizeof(Matrix));
    C->dimension = newDim;
    C->startRow = original->startRow + newDim;
    C->startColumn = original->startColumn;
    C->values = original->values;

    Matrix* D = (Matrix*) malloc(sizeof(Matrix));
    D->dimension = newDim;
    D->startRow = original->startRow + newDim;
    D->startColumn = original->startColumn + newDim;
    D->values = original->values;

    Matrix** res = (Matrix**) malloc(4 * sizeof(Matrix*));
    res[0] = A;
    res[1] = B;
    res[2] = C;
    res[3] = D;

    return res;
};

void combine(Matrix* Product, Matrix* TL, Matrix* TR, Matrix* BL, Matrix* BR){
    int startRow = Product->startRow;
    int startCol = Product->startColumn;
    int halfDim = TL->dimension;

    for (int i = 0; i < halfDim; i++){
        for (int j = 0; j < halfDim; j++){
            Product->values[i][j] = TL->values[i][j];
            
        }
    }

    for (int i = 0; i < halfDim; i++){
        for (int j = 0; j < halfDim; j++){
            Product->values[i][j + halfDim] = TR->values[i][j];
        }
    }

    for (int i = 0; i < halfDim; i++){
        for (int j = 0; j < halfDim; j++){
            Product->values[i + halfDim][j] = BL->values[i][j];
        }
    }

    for (int i = 0; i < halfDim; i++){
        for (int j = 0; j < halfDim; j++){
            Product->values[i + halfDim][j + halfDim] = BR->values[i][j];
        }
    }
}

Matrix* strassen(Matrix* m1, Matrix* m2){
    if (m1->dimension == N_0) {
        return conventionalMult(m1, m1);
    }
    
    Matrix** matrices1 = splitMatrices(m1);
    Matrix** matrices2 = splitMatrices(m2);
    Matrix* A = matrices1[0];
    Matrix* B = matrices1[1];
    Matrix* C = matrices1[2];
    Matrix* D = matrices1[3];
    Matrix* E = matrices2[0];
    Matrix* F = matrices2[1];
    Matrix* G = matrices2[2];
    Matrix* H = matrices2[3];

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

    /*
    AE +BG = −P2 +P4 +P5 +P6
    AF +BH = P1 +P2
    CE +DG = P3 +P4
    CF +DH = P1 −P3 +P5 +P7
    */

    Matrix* topLeft = addMatrices(addMatrices(P4, P2, true), addMatrices(P5, P6, false), false);
    Matrix* topRight = addMatrices(P1, P2, false);
    Matrix* bottomLeft = addMatrices(P3, P4, false);
    Matrix* bottomRight = addMatrices(addMatrices(P1, P3, true), addMatrices(P5, P7, false), false);

    Matrix* Product = initMatrix(m1->dimension);
    combine(Product, topLeft, topRight, bottomLeft, bottomRight);
    return Product;
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

void freeMatrix(Matrix* mat) {
	for (int i = 0; i < mat->dimension; i++) {
		free(mat->values[i]);
	}
	free(mat->values);
	free(mat);
}
