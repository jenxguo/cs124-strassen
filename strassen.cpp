#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
using namespace std;

int** strassen(int** m1, int** m2, int n);
int** bruteForce(int** m1, int** m2, int n);

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
        printf("Usage: ./randmst 0 dimension inputfile\n");
        return 1;
    }
    return 0;
}


