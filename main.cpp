#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <valarray>
#include <chrono>
#include "./model/LILMatrix.h"

using namespace std;
using namespace std::chrono;

#define epsilon 0.0001

void normalizeVector(vector<double> &vector) {
    double sum = 0;
    for (double i : vector)
        sum += i;
    for (double& i : vector)
        i = i / sum;
}

int main(int argc, char **argv) {

    if (argc != 3) {
       printf("Cantidad de parámetros incorrecta.\n");
       return 1;
    }



    // INITIALIZE
    char *inputFile = argv[1];
    double p = atof(argv[2]);

    ifstream fileInput;
    fileInput.open(inputFile);

    int pages;
    int links;
    fileInput >> pages >> links;



    // TIME START
    auto start = high_resolution_clock::now();

    cout << "____________________________________________________________________________" << endl;
    cout << "Definiendo matriz W..." << endl;

    LILMatrix W(pages, links,fileInput);

    cout << "Matriz W definida con dimensiones: " << W.getRows() << "x" << W.getCols() << endl << endl;


    cout << "Definiendo matriz D...";

    vector<double> D(pages);
    for (int i = 0; i < pages; i++) {
       double pageGrade = W.getPageGrade(i);

       if (pageGrade == 0) {
           D[i] = 0;
           continue;
       }
       double newValue = 1/pageGrade;
       D[i] = abs(newValue) >= epsilon ? newValue : 0;
    }

    cout << "Matriz D definida." << endl << endl;


    cout << "Definiendo vector e...";
    vector<double> e(pages, 1);

    cout << "Vector e definido." << endl;

    cout << "--------------------------------------------------" << endl;

    cout << "Iniciando page rank algorithm..." << endl << endl;

    // PAGE RANK
    W.multiplicationByScalar(p);


    W.multiplicationByDiagonalMatrix(D);

    W.identitySubtractSelf();
    W.gaussianElimination(e);
    normalizeVector(e);

    cout << endl;
    cout << "Terminando page rank algorithm...";

    // TIME END
    auto stop = high_resolution_clock::now();

    auto micro = duration_cast<microseconds>(stop - start);
    auto ms = duration_cast<milliseconds>(stop - start);
    auto s = duration_cast<seconds>(stop - start);
    auto min = duration_cast<minutes>(stop - start);



    // OUTPUT
    ofstream Output;
    Output.open(string(inputFile) + ".own.out");
    Output << p << endl;

    for (double value : e)
       Output << value << endl;

    cout << endl;
    cout << "--------------------------------------------------" << endl;
    cout << "INPUT FILE NAME:   " << string(inputFile) << endl;
    cout << "OUTPUT FILE NAME:   " << string(inputFile) + ".own.out" << endl;
    cout << "P VALUE:   " << p << endl;
    cout << "RESULT: " << e[0] << endl;
    double sum = e[0];
    for (int i = 1; i < e.size(); i++) {
        cout << "        " << e[i] << endl;
        sum += e[i];
    }
    cout << "SUM:    " << sum << endl;

    cout << "TIMER:   " << micro.count() << "us   |   " << ms.count() << "ms   |   " << s.count() << "s   |   " << min.count() << "min" << endl;
    cout << "____________________________________________________________________________" << endl;

    return 0;
}