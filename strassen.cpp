#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[]) {

    //checks that there are command line args
    if (argc != 4)
    {
        printf("Usage: ./randmst 0 dimension inputfile\n");
        return 1;
    }
    return 0;
}