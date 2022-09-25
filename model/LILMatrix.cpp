#include "LILMatrix.h"
#include <fstream>
#include <cassert>

/* CONSTRUCTORS */

LILMatrix::LILMatrix(int rows, int columns) : rows(rows), cols(columns), nz_elems(0) {}

LILMatrix::LILMatrix(int pages, int links, ifstream &file) {
    this->rows = pages;
    this->cols = pages;
    this->nz_elems = 0;
    for (int i = 0; i < links; i++) {
        int srcPage, dstPage;
        file >> srcPage >> dstPage;
        setValue(srcPage-1, dstPage-1, 1);
    }
}

LILMatrix::LILMatrix(const vector<vector<double>> &vectorMatrix) {
    this->rows = (int) (vectorMatrix.size());
    this->cols = (int) vectorMatrix[0].size();
    this->nz_elems = 0;

    for (int i = 0; i < rows; i++) {
        assert(vectorMatrix[i].size() == cols);
        for (int j = 0; j < cols; j++)
            setValue(i, j, vectorMatrix[i][j]);
    }
}



/* ACCESSING */

matrix_t::iterator LILMatrix::findRow(int target) {
    //for (auto row : matrix)
    //    if (row.first == target) return row;


    auto m_it = matrix.begin();
    while(m_it != matrix.end() && m_it->first < target) m_it++;
    return m_it;

}
row_elements_t::iterator LILMatrix::findColumn(int target, row_elements_t* row) {
    auto row_it = (*row).begin();
    while(row_it != (*row).end() && row_it->first < target) row_it++;
    return row_it;
}
// TODO la idea es sacar esta función en la entrega final porque empeora la performance
// Esto es solo para poder correr los test y verificar que no hay vectores con elementos de más
// Si el nz_elems se maneja con incrementos (++) y decrementos (--),
// no se sabe realmente cuántos elementos no nulos hay en la matriz
// Por lo tanto esta función los cuenta y reasigna el valor que nz_elems debería tener realmente.
int LILMatrix::countNzElems() {
    int ggwp = 0;
    for (auto& i : matrix)
        for (auto& j : i.second)
            ggwp++;
    return ggwp;
}

void LILMatrix::setValue(int targetRow, int targetColumn, double targetValue) {
    assert(targetRow >= 0 && targetRow <= rows - 1);
    assert(targetColumn >= 0 && targetColumn <= cols - 1);

    if(matrix.empty()) {
        if(targetValue == 0) return;
        matrix.push_back(row_t(targetRow, row_elements_t(0)));
        auto &row_elem = matrix.back();
        row_elem.second.emplace_back(targetColumn, targetValue);
        //nz_elems = countNzElems();
        nz_elems++;
        return;
    }

    value_t newElement = make_pair(targetColumn, targetValue); // elemento
    row_elements_t newRow_elem{newElement};
    row_t newMatrix_elem = make_pair(targetRow, newRow_elem); // fila

    double actualValue = getValue(targetRow, targetColumn);

    if (targetValue == 0) {
        if(actualValue == 0) return;

        auto m_it = findRow(targetRow);
        auto row_it = findColumn(targetColumn, &m_it->second);

        m_it->second.erase(row_it);

        //nz_elems = countNzElems();
        nz_elems--;
        return;
    }

    if (actualValue == 0) {
        auto m_it = findRow(targetRow);
        if(m_it->first != targetRow) matrix.insert(m_it, newMatrix_elem);
        else m_it->second.insert(findColumn(targetColumn, &m_it->second), newElement);

        //nz_elems = countNzElems();
        nz_elems++;
    } else {
        for(auto & i : matrix) { // sujeto a ver optimización
            if(i.first != targetRow) continue;
            for(auto & j : i.second) {
                if(j.first != targetColumn) continue;
                j.second = targetValue;
                return;
            }
        }
    }
}

double LILMatrix::getValue(int targetRow, int targetColumn) {
    assert(targetRow >= 0 && targetRow <= rows - 1);
    assert(targetColumn >= 0 && targetColumn <= cols - 1);

    auto m_it = findRow(targetRow);
    if(m_it == matrix.end() || targetRow != m_it->first) return 0;
    auto row_it = findColumn(targetColumn, &m_it->second);
    if(row_it == m_it->second.end() || targetColumn != row_it->first) return 0;

    return row_it->second;
}

double LILMatrix::getPageGrade(int targetRow) {
    auto m_it = findRow(targetRow);
    if (m_it == matrix.end() || targetRow != m_it->first) return 0;
    return (double) m_it->second.size();
    /*
    double res = 0;
    for (int i = 0; i < cols; i++) res += getValue(targetRow, i);
    return res;
    */
}

void LILMatrix::multiplicationByScalar(double scalar) {
    if (scalar == 0) {
        matrix = matrix_t(0);
        nz_elems = 0;
    }
    for (auto &i: matrix)
        for (auto &j: i.second)
            j.second = j.second * scalar;
}

void LILMatrix::multiplicationByDiagonalMatrix(vector<double> triangularMatrix) {
    for (auto &i : matrix)
        for (auto &j : i.second)
            setValue(i.first, j.first, triangularMatrix[j.first] * getValue(i.first, j.first));
}

void LILMatrix::identitySubtractSelf() {
    for (auto &i: matrix) {
        for (auto &j: i.second) {
            int row = i.first;
            int col = j.first;
            if (row == col) continue;
            setValue(row, col, (-1) * getValue(row, col));
        }
    }
    for (int k = 0; k < rows; k++) // TODO verificar si se puede meter en el ciclo anterior
        //setValue(k,k,1); // considering always zero values in the self diagonal //TODO a ver si esta optimización rompe
        setValue(k, k, 1 - getValue(k, k)); // considering nonzero values in the self diagonal}
}





void LILMatrix::gaussianElimination(vector<double> &b) {

    for (int di = 0; di < matrix.size(); di++) {
        double pivot = getValue(di, di);

        for (int ri = di + 1; ri < matrix.size(); ri++) {
            double rv = getValue(ri,di);
            if (rv == 0) continue;
            double factor = rv / pivot;

            for (int dei = 0; dei < matrix[di].second.size(); dei++) {

                int colCoord = matrix[di].second[dei].first;
                double pi = matrix[di].second[dei].second;

                rv = getValue(ri, colCoord);
                double newValue = rv - factor*pi;

                setValue(ri, colCoord, newValue);
            }
            b[ri] -= factor*b[di];
        }
    }

    for (int eci = (int) matrix.size() - 1; eci >= 0; eci--) {

        double sum = 0;
        for (int ti = 1; ti < matrix[eci].second.size(); ti++)
            sum += matrix[eci].second[ti].second * b[eci+ti];

        b[eci] = (b[eci] - sum) / matrix[eci].second[0].second;
    }
}





void LILMatrix::debug_structural_matrix() {
    cout << "-----" << endl;
    for (auto i : matrix) {
        for (auto j : i.second) {
            cout << "[" << i.first + 1 << "," << j.first + 1 << ": " << j.second << "]";
        }
        cout << endl;
    }
    cout << "-----" << endl;
}

void LILMatrix::debug_abstract_matrix() {
    cout << "-----" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << getValue(i,j) << " ";
        }
    }
    cout << "-----" << endl;
}

void LILMatrix::debug_ec_system() {
    cout << "-----" << endl;
    for (int i = 0; i < matrix.size(); i++) {
        cout << "equation " << i + 1 << ":   ";
        for (int j = 0; j < matrix[i].second.size(); j++) {
            cout << "(" << matrix[i].second[j].second << ") *x" << matrix[i].second[j].first+1 << "   ";
        }
        cout << endl;
    }
    cout << "-----" << endl;
}

int LILMatrix::getRows() const { return rows; }
int LILMatrix::getCols() const { return cols; }
int LILMatrix::nzElems() const { return nz_elems; }