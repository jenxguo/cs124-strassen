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

int numOpsForStrassen(int n0, int d);
int calcPadding(int d, int n_0);
Matrix* strassen(Matrix* m1, Matrix* m2, int flag, int n_0, Matrix** tempMatrices);
Matrix* conventionalMult(Matrix* m1, Matrix* m2);
int countTriangles(double p, int flag);

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
        //printf("Padding: %i\n", pad);

        A = initMatrix(pad);
        B = initMatrix(pad);
    }

    if (flag == 0 || flag == 1) {
        populateMatrices(A, B, d, inputfile);
    } else {
        generateRandomMatrix(A, d, 0);
        generateRandomMatrix(B, d, -1);
    }

    // printMat(A);
    // printf("\n");
    // printMat(B);
    // printf("\n");
    printf("starting strassen's\n");
    Matrix* C = strassen(A, B, flag, GLOBAL_N_0, tempMatrices);
    // printf("dimensions of strassen %i\n", C->dimension);
    if (flag != 1 && flag != 3) {
        C->dimension = d;
    }

    //printf("\n");
    //printMat(C);

    // printf("THE CORRECT ONE\n");
    Matrix* D = conventionalMult(A, B);
    // printf("dimensions of orrect %i\n", D->dimension);
    // printf("\n");
    if (flag != 1 && flag != 3) {
        D->dimension = d;
    }
    //printMat(D);


    /* Required Output for Autograder*/
    for (int v = 0; v < d; v++){
        printf("%i\n", C->values[v][v]);
    }

    if (checkCorrectness(C, D)) {
        printf("CORRECT\n");
    } else {
        printf("INCORRECT :( FUUUUU\n");
    }

    for (int i = 0; i < 18; i++){
        freeMatrix(tempMatrices[i]);
    }

    printf("%i\n", numOpsForStrassen(40, 233));

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
    myfile << "n_0, dimension, time1, time2, time3, time4, time5, average" << endl;

    for (int dim = 900; dim <= 900; dim += 51){
        Matrix** tempMatrices = (Matrix**) malloc(18 * sizeof(Matrix*));
        for (int i = 0; i < 18; i++){
            tempMatrices[i] = initMatrix(dim + 2);
        }

        for (int test_n_0 = 35; test_n_0 <= 85; test_n_0 += 4){
            myfile << test_n_0 << ", " << dim;

            // Test out strassen's algorithm 5 times
            duration<double, std::milli> totalTime;
            for (int i = 0; i < 5; i++){
                // Randomly create matrices A and B
                // int pad = calcPadding(dim, test_n_0);
                Matrix* A = initMatrix(dim);
                Matrix* B = initMatrix(dim);

                generateRandomMatrix(A, dim, 0);
                generateRandomMatrix(B, dim, 0);

                std::chrono::high_resolution_clock::time_point start = high_resolution_clock::now();
                Matrix* strassenRes = strassen(A, B, 1, test_n_0, tempMatrices);
                std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();
                duration<double, std::milli> ms_double = end - start;
                myfile << ", " << ms_double.count();
                if (i == 0){
                    totalTime = ms_double;
                } else {
                    totalTime += ms_double;
                }
                Matrix* convRes = conventionalMult(A, B);
                if (checkCorrectness(strassenRes, convRes)) {
                    printf("yeet\n");
                    // printMat(strassenRes);
                    // printf("\n");
                    // printMat(convRes);
                } else {
                    printf("try again it's gon be all right\n");
                    // printMat(strassenRes);
                    // printf("\n");
                    // printMat(convRes);
                }
                freeMatrix(A);
                freeMatrix(B);
                freeMatrix(strassenRes);
                // freeMatrix(convRes);
            }
            myfile << ", " << (totalTime / 5).count() << endl;
        }
        for (int i = 0; i < 18; i++){
            freeMatrix(tempMatrices[i]);
        }
    }
    
    myfile.close();
};

// main Strassen's algorithm to multiply matrices
Matrix* strassen(Matrix* m1, Matrix* m2, int flag, int n_0, Matrix** tempMatrices){
    
    int d = m1->dimension;
    if (d <= n_0) {
        //printf("\n new matrix sizes are %i and %i\n", m1->dimension, m2->dimension);
        return conventionalMult(m1, m2);
    }

    Matrix* newm1;
    Matrix* newm2;
    bool padded = false;

    // pad if odd d
    if ((flag == 1 || flag == 3 || flag == 4) && (d % 2) == 1) {
        // sizeUp(m1);
        // sizeUp(m2);
        newm1 = initMatrix(d+1);
        newm2 = initMatrix(d+1);
        copyMatrix(m1, newm1);
        copyMatrix(m2, newm2);
        d += 1;
        padded = true;
        // freeMatrix(m1);
        // freeMatrix(m2);
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

    
    int newDim = A->dimension; // newm1->dimension / 2;
    tempMatrices[0]->dimension = newDim;
    tempMatrices[0]->startRow = newDim;
    tempMatrices[0]->startColumn = newDim;
    tempMatrices[0] = addMatrices(F, H, tempMatrices[0], true);
    Matrix* P1 = strassen(A, tempMatrices[0], flag, n_0, tempMatrices);

    tempMatrices[1]->dimension = newDim;
    tempMatrices[1]->startRow = newDim;
    tempMatrices[1]->startColumn = newDim;
    tempMatrices[1] = addMatrices(A, B, tempMatrices[1], false);
    Matrix* P2 = strassen(tempMatrices[1], H, flag, n_0, tempMatrices);
    
    tempMatrices[2]->dimension = newDim;
    tempMatrices[2]->startRow = newDim;
    tempMatrices[2]->startColumn = newDim;
    tempMatrices[2] = addMatrices(C, D, tempMatrices[2], false);
    Matrix* P3 = strassen(tempMatrices[2], E, flag, n_0, tempMatrices);
       
    tempMatrices[3]->dimension = newDim;
    tempMatrices[3]->startRow = newDim;
    tempMatrices[3]->startColumn = newDim;
    Matrix* P4 = strassen(D, addMatrices(G, E, tempMatrices[3], true), flag, n_0, tempMatrices);

    tempMatrices[4]->dimension = newDim;
    tempMatrices[4]->startRow = newDim;
    tempMatrices[4]->startColumn = newDim;
    tempMatrices[5]->dimension = newDim;
    tempMatrices[5]->startRow = newDim;
    tempMatrices[5]->startColumn = newDim;
    tempMatrices[4] = addMatrices(A, D, tempMatrices[4], false);
    tempMatrices[5] = addMatrices(E, H, tempMatrices[5], false);
    Matrix* P5 = strassen(tempMatrices[4], tempMatrices[5], flag, n_0, tempMatrices);

    tempMatrices[6]->dimension = newDim;
    tempMatrices[6]->startRow = newDim;
    tempMatrices[6]->startColumn = newDim;
    tempMatrices[7]->dimension = newDim;
    tempMatrices[7]->startRow = newDim;
    tempMatrices[7]->startColumn = newDim;
    tempMatrices[6] = addMatrices(B, D, tempMatrices[6], true);
    tempMatrices[7] = addMatrices(G, H, tempMatrices[7], false);
    Matrix* P6 = strassen(tempMatrices[6], tempMatrices[7], flag, n_0, tempMatrices);

    tempMatrices[8]->dimension = newDim;
    tempMatrices[8]->startRow = newDim;
    tempMatrices[8]->startColumn = newDim;
    tempMatrices[9]->dimension = newDim;
    tempMatrices[9]->startRow = newDim;
    tempMatrices[9]->startColumn = newDim;
    tempMatrices[8] = addMatrices(C, A, tempMatrices[8], true);
    tempMatrices[9] = addMatrices(E, F, tempMatrices[9], false);
    Matrix* P7 = strassen(tempMatrices[8], tempMatrices[9], flag, n_0, tempMatrices);


    // free array holding split matrices
    free(matrices1);
    free(matrices2);

    Matrix* topLeft = addMatrices(addMatrices(P4, P2, tempMatrices[10], true), addMatrices(P5, P6, tempMatrices[11], false), tempMatrices[12], false);
    Matrix* topRight = addMatrices(P1, P2, tempMatrices[13], false);
    Matrix* bottomLeft = addMatrices(P3, P4, tempMatrices[14], false);
    Matrix* bottomRight = addMatrices(addMatrices(P1, P3, tempMatrices[15], true), addMatrices(P5, P7, tempMatrices[16], false), tempMatrices[17], false);

    Matrix* Product = initMatrix(d); 
    combine(Product, topLeft, topRight, bottomLeft, bottomRight);
    
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

// Combines 4 matrices to one big matrix
void combine(Matrix* Product, Matrix* TL, Matrix* TR, Matrix* BL, Matrix* BR){
    int halfDim = Product->dimension / 2;

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

// int countTriangles(double p, int flag) {
//     srand (static_cast <unsigned> (time(0)));
//     int size = 1024;
//     Matrix* graph = initMatrix(size);
//     for (int i = 0; i < size; i++) {
//         for (int j = 0; j < size; j++) {
//             int el;
//             if (rand() <  p * ((double)RAND_MAX + 1.0)) {
//                 el = 1;
//             } else {
//                 el = 0;
//             }
//             graph->values[i][j] = el;
//         }
//     }
//     Matrix* graphSquared = strassen(graph, graph, flag, GLOBAL_N_0);
//     Matrix* graphCubed = strassen(graphSquared, graph, flag, GLOBAL_N_0);
//     int diagonalSum = 0;
//     for (int i = 0; i < size; i++) {
//         diagonalSum += graphCubed->values[i][i];
//     }
//     return diagonalSum / 6;
// }

int numOpsForStrassen(int n0, int d) {
    if (d <= 1) {
        return 1;
    }
    if (d <= n0) {
        return (2 * pow(d, 3)) - pow(d, 2);
    }
    if (d % 2 == 1) {
        d+=1;
    }
    return 7*numOpsForStrassen(n0, d/2)+ 18*pow(d/2, 2);
}

// Frees matrix
void freeMatrix(Matrix* mat) {
	for (int i = 0; i < mat->dimension; i++) {
		free(mat->values[i]);
	}
	free(mat->values);
	free(mat);
}
