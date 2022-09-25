#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <valarray>
#include <chrono>
#include "./model/LILMatrix.h"

using namespace std;
using namespace std::chrono;


void normalizeVector(vector<double> &vector) {
    double sum = 0;
    for (double i : vector)
        sum += i;
    for (double& i : vector)
        i = i / sum;
}

void normalizeTester() {
    double numerador = -25;
    double denominador = 89;
    double frac = numerador / denominador;
    vector<double> arr = {frac, frac, frac, frac, 1};
    normalizeVector(arr);
    for (auto i : arr)
        cout << "Valor: " << i << endl;

    double sum = 0;
    for (double& i : arr)
        sum += i;
    cout << sum << endl;
}

int main(int argc, char **argv) {

    if (argc != 3) {
       printf("Cantidad de parámetros incorrecta.\n");
       return 1;
    }



    // TIME START
    auto start = high_resolution_clock::now();



    // INITIALIZE
    char *inputFile = argv[1];
    double p = atof(argv[2]);

    ifstream fileInput;
    fileInput.open(inputFile);

    int pages;
    int links;
    fileInput >> pages >> links;

    cout << "____________________________________________________________________________" << endl;
    cout << "Definiendo matriz W..." << endl;

    LILMatrix W(pages, links,fileInput);

    cout << "Matriz W definida con dimensiones: " << W.getRows() << "x" << W.getCols() << endl << endl;


    cout << "Definiendo matriz D..." << endl;

    vector<double> D(pages);
    for (int i = 0; i < pages; i++) {
       double pageGrade = W.getPageGrade(i);
       D[i] = pageGrade == 0 ? 0 : pageGrade;
    }

    cout << "Matriz D definida." << endl << endl;


    cout << "Definiendo vector e..." << endl;
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
    cout << "FILE NAME:   " << string(inputFile) << endl;
    cout << "P VALUE:   " << p << endl;
    cout << "RESULT: " << e[0] << endl;
    double sum = e[0];
    for (int i = 1; i < e.size(); i++) {
        cout << "        " << e[i] << endl;
        sum += e[i];
    }
    cout << "SUM:    " << sum << endl;

    cout << "TIMER:   " << micro.count() << "us   |   " << ms.count() << "ms   |   " << s.count() << "s   |   " << min.count() << "min" << endl;
    cout << "--------------------------------------------------" << endl;

    //normalizeTester();

    /*
    double numerador = -25;
    double denominador = 89;
    double frac = numerador / denominador;
    vector<double> arr = {frac, frac, frac, frac, 1};

    normalizeVector(arr);

    for (auto i : arr)
        cout << "Valor: " << i << endl;

    normalizeVector(arr);

    for (auto i : arr)
        cout << "Valor 2: " << i << endl;
    */

    return 0;
}














/* LEGACY
vector<double> z(pages);
for (int i = 0; i < pages; i++) {
   double pageGrade = W.getPageGrade(i);
   z[i] = pageGrade == 0 ? 1.0 / pages : (1 - p) / pages;
}
*/










/*
vector<int> arr = {1,2,3,4,5,6,7,8,9,10};
auto it = arr.begin();

for (int i = 0; i < arr.size(); i++) {
    if (arr[i] == 5) {
        arr.erase(it);
        i--;
    }
    arr[i]++;
}
for (auto& i : arr) {
    if (i == 5) arr.erase(it);
    i++;
}
for (auto& i : arr)
    cout << i << endl;
*/