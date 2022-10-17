#ifndef C_LILMATRIX_H
#define C_LILMATRIX_H

#include <iostream>
#include <vector>

using namespace std;

using ROW_NUMBER = int;
using COLUMN_NUMBER = int;
using ELEMENT_VALUE = double;

using value_t       = pair<COLUMN_NUMBER, ELEMENT_VALUE>;
using row_elements_t= vector<value_t>;
using row_t         = pair<ROW_NUMBER, row_elements_t>;
using matrix_t      = vector<row_t>;

class LILMatrix {
private:
    matrix_t matrix;
    int rows;
    int cols;
    int nz_elems; // the number of nonzero elements
    matrix_t::iterator findRow(int target);
    row_elements_t::iterator findColumn(int target, row_elements_t* row);

public:
    // Constructors
    LILMatrix(int rows, int columns);
    explicit LILMatrix(int pages, int links, ifstream &file);
    explicit LILMatrix(const vector<vector<double>> &vectorMatrix);

    // Accessing
    int getRows() const;
    int getCols() const;
    int nzElems() const;
    double getValue(int targetRow, int targetColumn);
    void setValue(int row, int column, double targetValue);

    // Actions
    double getPageGrade(int page);
    void multiplicationByScalar(double scalar);
    void multiplicationByDiagonalMatrix(vector<double> triangularMatrix);
    void identitySubtractSelf();
    void gaussianElimination(vector<double> &b);

    // Extra
    void debug_ec_system(const vector<double>& b);
    void debug_abstract_matrix();
    void debug_structural_matrix();
};


#endif
