#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>
using namespace std;

struct Matrix{
    int startRow;
    int startColumn;
    int dimension;
    int** values;
};

int calcPadding(int d, int n_0);
Matrix* strassen(Matrix* m1, Matrix* m2, int flag, int n_0);
Matrix* conventionalMult(Matrix* m1, Matrix* m2);
int countTriangles(double p, int flag);
int triangleTrial(int trials, double p, int flag);

// Matrix Functions
void generateRandomMatrix(Matrix* mat, int d, int offset);
Matrix* initMatrix(int d);
void copyMatrix(Matrix* oldMat, Matrix* newMat);
Matrix** splitMatrices(Matrix* original);
void combine(Matrix* Product, Matrix* TL, Matrix* TR, Matrix* BL, Matrix* BR);
Matrix* addMatrices(Matrix* A, Matrix* B, Matrix* res, bool subtract);
void populateMatrices(Matrix* A, Matrix* B, int d, char* inputfile);
void printMat(Matrix* mat);
void freeMatrix(Matrix* mat);
bool checkCorrectness(Matrix* A, Matrix* B);
void sizeUp(Matrix* m);

void runTests();

// Transition value of n to start using conventional mult alg
int GLOBAL_N_0 = 50;

int main(int argc, char *argv[]) {

    //checks that there are command line args
    if (argc != 4)
    {
        cout << "Usage: ./strassen 0 dimension inputfile\n";
        return 1;
    }

    // FLAG CODE
    // 0: input file, initial padding method
    // 1: input file, recursive padding method
    // 2: generate random matrices, initial padding method
    // 3: generate random matrices, recursive padding method
    // 4: just testing!
    int flag = (int) strtol(argv[1], NULL, 10);
    int d = (int) strtol(argv[2], NULL, 10);
    char *inputfile = argv[3];

    /* Just testing */
    if (flag == 4) {
        runTests();
        return 0;   
    }

    Matrix** tempMatrices = (Matrix**) malloc(18 * sizeof(Matrix*));
    for (int i = 0; i < 18; i++){
        tempMatrices[i] = initMatrix(d);
    }

    Matrix* A;
    Matrix* B;
    // Recursive Padding Method
    if (flag == 1 || flag == 3) {
        A = initMatrix(d);
        B = initMatrix(d);
    }
    // Initial Padding Method
    else {
        int pad = calcPadding(d, GLOBAL_N_0);
        A = initMatrix(pad);
        B = initMatrix(pad);
    }

    if (flag == 0 || flag == 1) {
        populateMatrices(A, B, d, inputfile);
    } else {
        generateRandomMatrix(A, d, 0);
        generateRandomMatrix(B, d, -1);
    }

    Matrix* C = strassen(A, B, flag, GLOBAL_N_0);
    if (flag != 1 && flag != 3) {
        C->dimension = d;
    }

    Matrix* D = conventionalMult(A, B);
    if (flag != 1 && flag != 3) {
        D->dimension = d;
    }

    /* Required Output for Autograder*/
    for (int v = 0; v < d; v++){
        printf("%i\n", C->values[v][v]);
    }

    return 0;
}

/* Find experimental value of N_0 */
void runTests() {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // Open output file
    ofstream myfile;
    myfile.open("data.txt", ofstream::app);
    myfile << "n_0, dimension, time1, time2, time3, time4, time5, time6, time7, time8, time9, time10, average" << endl;

    for (int dim = 1300; dim <= 1300; dim += 51){

        for (int test_n_0 = 115; test_n_0 <= 300; test_n_0 += 10){
            myfile << test_n_0 << ", " << dim;

            // Test for 10 trials
            duration<double, std::milli> totalTime;
            for (int i = 0; i < 10; i++){
                // Randomly create matrices A and B
                Matrix* A = initMatrix(dim);
                Matrix* B = initMatrix(dim);

                generateRandomMatrix(A, dim, 0);
                generateRandomMatrix(B, dim, 0);

                std::chrono::high_resolution_clock::time_point start = high_resolution_clock::now();
                Matrix* strassenRes = strassen(A, B, 1, test_n_0);
                std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();
                duration<double, std::milli> ms_double = end - start;
                myfile << ", " << ms_double.count();
                if (i == 0){
                    totalTime = ms_double;
                } else {
                    totalTime += ms_double;
                }
                // Matrix* convRes = conventionalMult(A, B);
                freeMatrix(A);
                freeMatrix(B);
                freeMatrix(strassenRes);
            }
            myfile << ", " << (totalTime / 10).count() << endl;
        }
    }
    
    myfile.close();
};

// main Strassen's algorithm to multiply matrices
Matrix* strassen(Matrix* m1, Matrix* m2, int flag, int n_0){
    int d = m1->dimension;
    if (d <= n_0) {
        return conventionalMult(m1, m2);
    }

    Matrix* newm1;
    Matrix* newm2;
    bool padded = false;

    // pad if odd d
    if ((flag == 1 || flag == 3 || flag == 4) && (d % 2) == 1) {
        newm1 = initMatrix(d+1);
        newm2 = initMatrix(d+1);
        copyMatrix(m1, newm1);
        copyMatrix(m2, newm2);
        d += 1;
        padded = true;
    } 
    else {
        newm1 = m1;
        newm2 = m2;
    }
    
    
    Matrix** matrices1 = splitMatrices(newm1);
    Matrix** matrices2 = splitMatrices(newm2);
    Matrix* A = matrices1[0];
    Matrix* B = matrices1[1];
    Matrix* C = matrices1[2];
    Matrix* D = matrices1[3];
    Matrix* E = matrices2[0];
    Matrix* F = matrices2[1];
    Matrix* G = matrices2[2];
    Matrix* H = matrices2[3];
    
    if (C->dimension != D->dimension) {
        printf("\nhelp\n");
        printf("newm1 is %i and newm2 is %i\n", newm1->dimension, newm2->dimension);
        printf("the dimensions are as follows: %i %i %i %i %i %i %i %i\n", A->dimension, B->dimension, C->dimension, D->dimension, E->dimension, F->dimension, G->dimension, H->dimension);
    }

    Matrix* temp1 = initMatrix(A->dimension);
    Matrix* temp2 = initMatrix(A->dimension);
    
    addMatrices(F, H, temp1, true);
    Matrix* P1 = strassen(A, temp1, flag, n_0);

    addMatrices(A, B, temp1, false);
    Matrix* P2 = strassen(temp1, H, flag, n_0);
    
    addMatrices(C, D, temp1, false);
    Matrix* P3 = strassen(temp1, E, flag, n_0);

    addMatrices(G, E, temp1, true);
    Matrix* P4 = strassen(D, temp1, flag, n_0);

    addMatrices(A, D, temp1, false);
    addMatrices(E, H, temp2, false);
    Matrix* P5 = strassen(temp1, temp2, flag, n_0);

    addMatrices(B, D, temp1, true);
    addMatrices(G, H, temp2, false);
    Matrix* P6 = strassen(temp1, temp2, flag, n_0);

    addMatrices(C, A, temp1, true);
    addMatrices(E, F, temp2, false);
    Matrix* P7 = strassen(temp1, temp2, flag, n_0);


    // free array holding split matrices
    free(matrices1);
    free(matrices2);

    // combine matrices by putting into prodct

    Matrix* Product = initMatrix(d);
    // top left
    Product->dimension = A->dimension;
    addMatrices(addMatrices(P4, P2, temp1, true), addMatrices(P5, P6, temp2, false), Product, false);

    // top right
    Product->startColumn = A->dimension;
    addMatrices(P1, P2, Product, false);

    // bottom right
    Product->startRow = A->dimension;
    addMatrices(addMatrices(P1, P3, temp1, true), addMatrices(P5, P7, temp2, false), Product, false);
    
    // bottom left
    Product->startColumn = 0;
    addMatrices(P3, P4, Product, false);

    Product->dimension = d;
    Product->startRow = 0;

    freeMatrix(temp1);
    freeMatrix(temp2);
    
    // get rid of extra 0 dimension if padded
    if ((flag == 1 || flag == 3 || flag == 4) && padded) {
        Product->dimension--;
    }
    
    return Product;
};


// Adds matrix A and B with flag for subtraction if subtract = true
Matrix* addMatrices(Matrix* A, Matrix* B, Matrix* res, bool subtract) {
    if (A->dimension != B->dimension) {
        printf("Error with adding matrices: %i vs %i\n", A->dimension, B->dimension);
        return NULL;
    } 

    int d = A->dimension;
    int ARidx = A->startRow;
    int ACidx = A->startColumn;
    int BRidx = B->startRow;
    int BCidx = B->startColumn;
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (subtract) {
                res->values[res->startRow + i][res->startColumn + j] = A->values[ARidx + i][ACidx + j] - B->values[BRidx + i][BCidx + j];
            } else {
                res->values[res->startRow + i][res->startColumn + j] = A->values[ARidx + i][ACidx + j] + B->values[BRidx + i][BCidx + j];
            }
        }
    }
    return res;
}

// Compares matrix A to B and returns true if identical
bool checkCorrectness(Matrix* A, Matrix* B) {
    int d = A->dimension;
    if (d != B->dimension) {
        return false;
    }
    int ARidx = A->startRow;
    int ACidx = A->startColumn;
    int BRidx = B->startRow;
    int BCidx = B->startColumn;
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (A->values[ARidx + i][ACidx + j] != B->values[BRidx + i][BCidx + j]) {
                return false;
            }
        }
    }
    return true;
}

// Calculate dimensions for padded matrix in initial padding method
int calcPadding(int d, int n_0) {
    if (d <= n_0) {
        return d;
    }
    int n = d;
    while (n > n_0) {
        n = (n / 2) + ((n % 2) != 0) ;
    }
    while (n < d) {
        n = n*2;
    }
    return n;
}

// Initialize empty matrix object of 0s with dimension d
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

// Copy values from old matrix to new matrix
void copyMatrix(Matrix* oldMat, Matrix* newMat) {
    int d = oldMat->dimension;
    if (newMat->dimension < d) {
        printf("Error copying matrices\n");
        return;
    }
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            newMat->values[i][j] = oldMat->values[oldMat->startRow + i][oldMat->startColumn + j];
        }
    }
    int newDim = newMat->dimension;
    for (int k = 0; k < newDim; k++){
        newMat->values[k][newDim-1] = 0;
        newMat->values[newDim-1][k] = 0;
    }
}

// Populate values into two matrices from input file
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

// Prints matrix to console
void printMat(Matrix* mat) {
    int d = mat->dimension;
    for (int i = mat->startRow; i < mat->startRow + d; i++) {
        for (int j = mat->startColumn; j < mat->startColumn + d; j++){
            printf("%i ", mat->values[i][j]);
        }
        printf("\n");
    }
}


// Returns array of 4 split matrices from original matrix
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

// Conventional algorithm to multiply matrices
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

// Generates random matrix of dimension d of values {0, 1, 2} or {-1, 0, 1} if offset = -1
void generateRandomMatrix(Matrix* mat, int d, int offset) {
    // Reseed
    srand (static_cast <unsigned> (time(0)));
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            int el = rand() % 3 + offset;
            mat->values[i][j] = el;
        }
    }
}

// generates random 1024 x 1024 matrix representation of graph and calculates number of triangles
int countTriangles(double p, int flag) {
    srand (static_cast <unsigned> (time(0)));
    int size = 1024;
    Matrix* graph = initMatrix(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int el;
            if (rand() <  p * ((double)RAND_MAX + 1.0)) {
                el = 1;
            } else {
                el = 0;
            }
            graph->values[i][j] = el;
        }
    }
    Matrix* graphSquared = strassen(graph, graph, flag, GLOBAL_N_0);
    Matrix* graphCubed = strassen(graphSquared, graph, flag, GLOBAL_N_0);
    int diagonalSum = 0;
    for (int i = 0; i < size; i++) {
        diagonalSum += graphCubed->values[i][i];
    }
    return diagonalSum / 6;
}

// Finds average number of triangles for certain number of trials
int triangleTrial(int trials, double p, int flag) {
    int sum = 0;
    for (int i = 0; i < trials; i++) {
        sum += countTriangles(p, flag);
    }
    return sum / trials;
}

// Frees matrix
void freeMatrix(Matrix* mat) {
	for (int i = 0; i < mat->dimension; i++) {
		free(mat->values[i]);
	}
	free(mat->values);
	free(mat);
}
