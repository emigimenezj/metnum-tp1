#include "LILMatrix.h"
#include <fstream>
#include <cassert>

#define epsilon 0.0001


/* CONSTRUCTORS */
LILMatrix::LILMatrix(int rows, int columns) : rows(rows), cols(columns), nz_elems(0) {}

LILMatrix::LILMatrix(int pages, int links, ifstream &file) {
    this->rows = pages;
    this->cols = pages;
    this->nz_elems = 0;
    for (int i = 0; i < links; i++) {
        int srcPage, dstPage;
        file >> srcPage >> dstPage;
        setValue(dstPage-1, srcPage-1, 1);
    }
}

LILMatrix::LILMatrix(const vector<vector<double>> &vectorMatrix) {
    this->rows = (int) vectorMatrix.size();
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
    auto m_it = matrix.begin();
    while(m_it != matrix.end() && m_it->first < target) m_it++;
    return m_it;
}
row_elements_t::iterator LILMatrix::findColumn(int target, row_elements_t* row) {
    auto row_it = (*row).begin();
    while (row_it != (*row).end() && row_it->first < target) row_it++;
    return row_it;
}

void LILMatrix::setValue(int targetRow, int targetColumn, double targetValue) {
    assert(targetRow >= 0 && targetRow <= rows - 1);
    assert(targetColumn >= 0 && targetColumn <= cols - 1);

    if (matrix.empty()) {
        if (targetValue == 0) return;
        matrix.push_back(row_t(targetRow, row_elements_t(0)));
        auto &row_elem = matrix.back();
        row_elem.second.emplace_back(targetColumn, targetValue);

        nz_elems++;
        return;
    }

    value_t newElement = make_pair(targetColumn, targetValue); // new row element: pair<col, value>
    row_elements_t newRow_elem{newElement};
    row_t newMatrix_elem = make_pair(targetRow, newRow_elem); // new row: pair<row, elements>

    double actualValue = getValue(targetRow, targetColumn);

    if (targetValue == 0) {
        if (actualValue == 0) return;

        auto m_it = findRow(targetRow);
        auto row_it = findColumn(targetColumn, &m_it->second);

        m_it->second.erase(row_it);

        nz_elems--;
        return;
    }

    if (actualValue == 0) {
        auto m_it = findRow(targetRow);

        if (m_it->first != targetRow)
            matrix.insert(m_it, newMatrix_elem);
        else
            m_it->second.insert(findColumn(targetColumn, &m_it->second), newElement);

        nz_elems++;
    } else {
        for(auto & i : matrix) {
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

double LILMatrix::getPageGrade(int targetColumn) {

    double res = 0;
    for (auto &i: matrix)
        for (auto &j: i.second)
            if (j.first == targetColumn) res += 1;

    return res;
}

void LILMatrix::multiplicationByScalar(double scalar) {
    if (scalar == 0) {
        matrix = matrix_t(0);
        nz_elems = 0;
    }
    for (auto &i: matrix)
        for (auto &j: i.second)
            j.second *= scalar;
}

void LILMatrix::multiplicationByDiagonalMatrix(vector<double> triangularMatrix) {
    for (auto &i : matrix)
        for (auto &j : i.second) {
            double newValue = j.second * triangularMatrix[j.first];
            setValue(i.first, j.first, abs(newValue) >= epsilon ? newValue : 0);
        }
}

void LILMatrix::identitySubtractSelf() {
    multiplicationByScalar(-1);
    for (int di = 0; di < rows; di++)
        setValue(di, di, 1);
}





void LILMatrix::gaussianElimination(vector<double> &b) {
    for (int di = 0; di < matrix.size(); di++) {
        for (int ri = di + 1; ri < matrix.size(); ri++) {
            if (matrix[ri].second[0].first != di) continue;

            double dv = matrix[di].second[0].second;
            double rv = matrix[ri].second[0].second;
            double factor = rv / dv;

            int rei = 0;
            for (int dei = 0; dei < matrix[di].second.size(); dei++) {

                auto de = matrix[di].second[dei];
                auto re = matrix[ri].second[rei];
                dv = matrix[di].second[dei].second;
                rv = matrix[ri].second[rei].second;

                double newValue = rv - factor * dv;

                if (de.first < re.first) { // Caso 1
                    newValue = - factor * dv;
                    if (abs(newValue) <= epsilon) continue;

                    value_t newElement = make_pair(de.first, newValue);
                    matrix[ri].second.insert(matrix[ri].second.begin() + rei, newElement);

                    rei++;
                    nz_elems++;
                    continue;
                }

                if (de.first == re.first) { // Caso 2

                    if (abs(newValue) <= epsilon) { // Caso 2.1
                        matrix[ri].second.erase(matrix[ri].second.begin()+rei);

                        nz_elems--;
                        continue;
                    }

                    // Caso 2.2
                    matrix[ri].second[rei].second = newValue;
                    continue;
                }

                // Caso 3
                rei++;
                if (rei >= matrix[ri].second.size()) { // Caso 3.1
                    while (dei < matrix[di].second.size()) {
                        newValue = - factor * dv;
                        if (abs(newValue) <= epsilon) {
                            dei++;
                        } else {
                            value_t newElement = make_pair(matrix[di].second[dei].first, newValue);
                            matrix[ri].second.push_back(newElement);

                            nz_elems++;
                            dei++;
                        }
                    }
                    break;
                }
                dei--;
            }
            b[ri] -= factor * b[di];
        }
    }

    for (int eci = (int) matrix.size() - 1; eci >= 0; eci--) {
        double sum = 0;
        for (int ti = 1; ti < matrix[eci].second.size(); ti++)
            sum += matrix[eci].second[ti].second * b[matrix[eci].second[ti].first];

        b[eci] =  (b[eci] - sum) / matrix[eci].second[0].second;
    }
}





void LILMatrix::debug_structural_matrix() {
    cout << "----- STRUCTURAL MATRIX" << endl;
    for (auto const& i : matrix) {
        for (auto j : i.second) {
            cout << "[" << i.first + 1 << "," << j.first + 1 << ": " << j.second << "]   ";
        }
        cout << endl;
    }
    cout << "-----" << endl;
}

void LILMatrix::debug_abstract_matrix() {
    cout << "----- ABSTRACT MATRIX" << endl;
    for (int i = 0; i < rows; i++) {
        cout << "| ";
        for (int j = 0; j < cols; j++) {
            cout << getValue(i,j) << " | ";
        }
        cout << endl;
    }
    cout << "-----" << endl;
}

void LILMatrix::debug_ec_system(const vector<double>& b) {
    cout << "----- EC SYSTEM" << endl;
    for (int i = 0; i < matrix.size(); i++) {
        cout << "equation " << i + 1 << ":   ";
        for (auto & j : matrix[i].second) {
            cout << "(" << j.second << ") *x" << j.first+1 << "   ";
        }
        cout << " = " << b[i] << endl;
    }
    cout << "-----" << endl;
}

int LILMatrix::getRows() const { return rows; }
int LILMatrix::getCols() const { return cols; }
int LILMatrix::nzElems() const { return nz_elems; }