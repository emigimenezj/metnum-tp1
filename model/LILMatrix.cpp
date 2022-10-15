#include "LILMatrix.h"
#include <fstream>
#include <cassert>

#define epsilon 0.0000000000000001


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
    /*
    cout << "INIT" << endl;
    debug_abstract_matrix();
    debug_structural_matrix();
     */
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

    /*
    int mid;
    int low = 0;
    int high = (int) matrix.size();

    while (low < high) {
        mid = (low + high) / 2;
        if (target > matrix[mid].first) low = mid + 1;
        else high = mid;
    }

    if (low < matrix.size() && matrix[low].first < target) low ++;

    return low;
    */
}
row_elements_t::iterator LILMatrix::findColumn(int target, row_elements_t* row) {

    auto row_it = (*row).begin();
    while(row_it != (*row).end() && row_it->first < target) row_it++;
    return row_it;


/*
    int mid;
    int low = 0;
    int high = (int) row.size();

    while (low < high) {
        mid = (low + high) / 2;
        if (target > row[mid].first) low = mid + 1;
        else high = mid;
    }

    if (low < row.size() && row[low].first < target) low ++;

    return low;
*/

    /*
    int mid;
    int low = 0;
    int high = (int) row.size();
    while (low <= high) {
        mid = (low + high) / 2;
        if (row[mid].first == target) return mid;
        if (row[mid].first < target) low = mid + 1;
        else high = mid - 1;
    }
    return -1;
    */

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

    if (matrix.empty()) {
        if (targetValue == 0) return;
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
        if (actualValue == 0) return;

        auto m_it = findRow(targetRow);
        auto row_it = findColumn(targetColumn, &m_it->second);

        m_it->second.erase(row_it);

        //nz_elems = countNzElems();
        nz_elems--;
        return;
    }

    if (actualValue == 0) {
        auto m_it = findRow(targetRow);
        if (m_it->first != targetRow) matrix.insert(m_it, newMatrix_elem);
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

double LILMatrix::getPageGrade(int targetColumn) {

    double res = 0;
    for (auto &i: matrix)
        for (auto &j: i.second)
            if (j.first == targetColumn) res += 1;

    return res;
    /*
    auto m_it = findRow(targetRow);
    if (m_it == matrix.end() || targetRow != m_it->first) return 0;
    return (double) m_it->second.size();*/

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
            j.second *= scalar;
    //cout << "multiplication by scalar" << endl;
    //debug_abstract_matrix();
    //debug_structural_matrix();
}

void LILMatrix::multiplicationByDiagonalMatrix(vector<double> triangularMatrix) {
    for (auto &i : matrix)
        for (auto &j : i.second) {
            double newValue = j.second * triangularMatrix[j.first];
            setValue(i.first, j.first, abs(newValue) >= epsilon ? newValue : 0);
        }
    //cout << "mult. by diagonal" << endl;
    //debug_abstract_matrix();
    //debug_structural_matrix();
}

void LILMatrix::identitySubtractSelf() {
    multiplicationByScalar(-1);
    for (int di = 0; di < rows; di++)
        setValue(di, di, 1);
}





void LILMatrix::gaussianElimination(vector<double> &b) {

    cout << "TRIANGULANDO: ";
/*
    for (int di = 0; di < matrix.size(); di++) {
        for (int ri = di+1; ri < matrix.size(); ri++) {
            if (matrix[ri].second[0].first != di) continue;

            auto actualDiag = matrix[di].second;
            auto actualRow = matrix[ri].second;
            auto actualRowIndex = matrix[ri].first;

            double factor = actualRow[0].second / actualDiag[0].second;

            int actualRowElementIndex = 0;
            for (auto& diagElem : actualDiag) {
                auto rowElem = matrix[ri].second[actualRowElementIndex];

                if (diagElem.first < rowElem.first) {
                    auto tmp = - factor*diagElem.second;
                    auto newValue = abs(tmp) < epsilon ? 0 : tmp;
                    setValue(actualRowIndex,diagElem.first, newValue);
                    actualRowElementIndex++;
                } else if (diagElem.first == rowElem.first) {
                    auto tmp = rowElem.second - factor * diagElem.second;
                    auto newValue = abs(tmp) < epsilon ? 0 : tmp;
                    setValue(actualRowIndex,diagElem.first,newValue);
                } else {
                    actualRowElementIndex++;
                    if (actualRowElementIndex >= actualRow.size()) {
                        auto tmp = - factor*diagElem.second;
                        auto newValue = abs(tmp) < epsilon ? 0 : tmp;
                        setValue(actualRowIndex, diagElem.first, newValue);
                    }
                }
            }
        }
    }
*/



    for (int di = 0; di < matrix.size(); di++) {
        cout << di << " ";
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

                if (de.first < re.first) {
                    newValue = - factor * dv;
                    if (abs(newValue) <= epsilon) continue;

                    //cout << di << "," << ri << "," << dei << "," << rei << "|" << re.first << endl;
                    //setValue(matrix[ri].first, re.first, newValue);
                    value_t newElement = make_pair(de.first, newValue);
                    matrix[ri].second.insert(matrix[ri].second.begin() + rei, newElement);

                    rei++;
                    nz_elems++;
                    continue;
                }

                if (de.first == re.first) {

                    if (abs(newValue) <= epsilon) { // agregamos || dei == 0
                        matrix[ri].second.erase(matrix[ri].second.begin()+rei);

                        nz_elems--;
                        continue;
                    }

                    matrix[ri].second[rei].second = newValue;
                    continue;
                }

                //while (rei < matrix[ri].second.size() && re.first < de.first) rei++;
                rei++;
                if (rei >= matrix[ri].second.size()) {
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
            //cout << endl;

//            double tmp = b[ri] - factor * b[di];
//            if (abs(tmp) < epsilon) b[ri] = 0;
//            else b[ri] = tmp;

            b[ri] -= factor * b[di];
        }
    }

    /*
    for (const auto& fila : matrix)
        if (fila.second.size() >= 400) cout << fila.first << "|" << fila.second.size() << endl;
    */

    /*
    matrix[500].second[0].first = 2;
    for (auto & i : matrix)
        if (i.first != i.second[0].first)
            cout << "[" << i.first << " " << i.second[0].first << "]" << endl;
    */

    //debug_structural_matrix();
    //debug_abstract_matrix();
    //debug_ec_system(b);

    /*cout << endl;
    for (int i = 1997; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].second.size(); j++) {
            cout << "[" << i << "]" << matrix[i].second[j].first << " ";
        }
        cout << b[i] << endl;
    }*/

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










/* LEGACY - gauss
                int colCoord = matrix[di].second[dei].first;
                double pv = matrix[di].second[dei].second;
                rv = getValue(ri, colCoord);
                double newValue = rv - factor * pv;
                cout << "[ (f)" << ri << "     (c)" << colCoord << "     (v) " << rv << "     (nv) " << newValue << " ]";
                cout << " ------------------------- " << di << endl;
                if (newValue == rv) continue;
                setValue(ri, colCoord, newValue);
*/



/* setValue with INT
    if (targetValue == 0) {
        if(actualValue == 0) return;

        auto row = matrix[findRow(targetRow)].second;
        auto rei = findColumn(targetColumn, row);

        row.erase(row.begin()+rei);

        //nz_elems = countNzElems();
        nz_elems--;
        return;
    }

    int ri = findRow(targetRow);

    if(ri != targetRow) {
        matrix.insert(matrix.begin() + ri, newMatrix_elem);
        //nz_elems = countNzElems();
        nz_elems++;
        return;
    }

    auto row = matrix[ri].second;
    int rei = findColumn(targetColumn, row);

    if (rei != targetColumn) {
        row.insert(row.begin()+rei, newElement);
        //nz_elems = countNzElems();
        nz_elems++;
        return;
    }

    row[rei].second = targetValue;*/


/* getValue with INT
    int ri = findRow(targetRow);
    if (ri != targetRow) return 0;
    int rei = findColumn(targetColumn, matrix[ri].second);
    if (rei != targetColumn) return 0;

    return matrix[ri].second[rei].first;*/

/* pageGrade with INT
    int ri = findRow(targetRow);
    if (ri != targetRow) return 0;
    return (double) matrix[ri].second.size();
    */









/* getValue with ITERATORS
    auto m_it = findRow(targetRow);
    if(m_it == matrix.end() || targetRow != m_it->first) return 0;
    auto row_it = findColumn(targetColumn, &m_it->second);
    if(row_it == m_it->second.end() || targetColumn != row_it->first) return 0;

    return row_it->second;
    */



/* setValue with ITERATORS
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
*/

/* pageGrade with ITERATORS
    auto m_it = findRow(targetRow);
    if (m_it == matrix.end() || targetRow != m_it->first) return 0;
    return (double) m_it->second.size();

    double res = 0;
    for (int i = 0; i < cols; i++) res += getValue(targetRow, i);
    return res;
    */