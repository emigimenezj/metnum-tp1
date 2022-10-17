#include <gtest/gtest.h>
#include "../model/LILMatrix.h"

/* CONSTRUCTORS */

TEST(Matrix_accessing, new_matrix_with_raw_constructor) {
    LILMatrix matrix(2,2);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_accessing, new_matrix_with_vector_constructor) {
    LILMatrix matrix({
        {1,2,3},
        {4,0,6},
        {7,8,9}
    });
    EXPECT_EQ(matrix.getRows(), 3);
    EXPECT_EQ(matrix.getCols(), 3);
    EXPECT_EQ(matrix.nzElems(), 8);
}





/* ACCESSING */

TEST(Matrix_accessing, empty_matrix_only_has_zeros) {
    LILMatrix matrix(2,2);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
}
TEST(Matrix_accessing, can_add_a_value) {
    LILMatrix matrix(2,2);
    matrix.setValue(1,1, 0.5);
    EXPECT_EQ(matrix.getValue(1,1), 0.5);
    EXPECT_EQ(matrix.nzElems(), 1);
}
TEST(Matrix_accessing, add_a_value_do_not_change_other_positions) {
    LILMatrix matrix(2,2);
    matrix.setValue(0,0, 0.5);
    EXPECT_EQ(matrix.getValue(0,0), 0.5);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
}
TEST(Matrix_accessing, can_remove_a_value) {
    LILMatrix matrix(2,2);
    matrix.setValue(0,0, 0.5);
    matrix.setValue(0,0, 0);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_accessing, remove_a_value_do_not_change_other_positions) {
    LILMatrix matrix(2,2);
    matrix.setValue(0,0, 0.5);
    matrix.setValue(0,0,0);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
}
TEST(Matrix_accessing, can_add_multiples_values) {
    LILMatrix matrix(5,5);
    matrix.setValue(0,0, 0.5);
    matrix.setValue(0,1, 1.5);
    matrix.setValue(0,2, 2.5);
    matrix.setValue(0,3, 3.5);
    matrix.setValue(0,4, 4.5);
    EXPECT_EQ(matrix.getValue(0,0), 0.5);
    EXPECT_EQ(matrix.getValue(0,1), 1.5);
    EXPECT_EQ(matrix.getValue(0,2), 2.5);
    EXPECT_EQ(matrix.getValue(0,3), 3.5);
    EXPECT_EQ(matrix.getValue(0,4), 4.5);
}
TEST(Matrix_accessing, can_remove_multiples_values) {
    LILMatrix matrix(5,5);
    matrix.setValue(0,0, 0.5);
    matrix.setValue(0,1, 1.5);
    matrix.setValue(0,2, 2.5);
    matrix.setValue(0,3, 3.5);
    matrix.setValue(0,4, 4.5);
    matrix.setValue(0,0, 0);
    matrix.setValue(0,1, 0);
    matrix.setValue(0,2, 0);
    matrix.setValue(0,3, 0);
    matrix.setValue(0,4, 0);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(0,2), 0);
    EXPECT_EQ(matrix.getValue(0,3), 0);
    EXPECT_EQ(matrix.getValue(0,4), 0);
}
TEST(Matrix_accessing, can_overwrite_a_value) {
    LILMatrix matrix(2,2);
    matrix.setValue(0,0, 0.5);
    matrix.setValue(0,0,99.5);
    EXPECT_EQ(matrix.getValue(0,0), 99.5);
    EXPECT_EQ(matrix.nzElems(), 1);
}
TEST(Matrix_accessing, remove_the_nothing_do_anything) {
    LILMatrix matrix(2,2);
    matrix.setValue(0,0, 0);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.nzElems(), 0);
}





/* ACTIONS */

TEST(Matrix_page_grade, page_grade_of_row_with_full_values) {
    LILMatrix matrix({{1,2},{3,4}});
    EXPECT_EQ(matrix.getPageGrade(0), 2);
    EXPECT_EQ(matrix.getPageGrade(1), 2);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 4);
}
TEST(Matrix_page_grade, page_grade_of_row_with_some_values) {
    LILMatrix matrix({{1,0},{0,4}});
    EXPECT_EQ(matrix.getPageGrade(0), 1);
    EXPECT_EQ(matrix.getPageGrade(1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_page_grade, page_grade_of_row_without_values) {
    LILMatrix matrix({{1,2},{0,0}});
    EXPECT_EQ(matrix.getPageGrade(0), 1);
    EXPECT_EQ(matrix.getPageGrade(1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_scalar_multiply, multiplicate_matrix_with_full_values) {
    LILMatrix matrix({{1,2},{3,4}});
    matrix.multiplicationByScalar(2);
    EXPECT_EQ(matrix.getValue(0,0), 2);
    EXPECT_EQ(matrix.getValue(0,1), 4);
    EXPECT_EQ(matrix.getValue(1,0), 6);
    EXPECT_EQ(matrix.getValue(1,1), 8);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 4);
}
TEST(Matrix_scalar_multiply, multiplicate_matrix_with_some_values) {
    LILMatrix matrix({{0,2},{3,0}});
    matrix.multiplicationByScalar(2);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 4);
    EXPECT_EQ(matrix.getValue(1,0), 6);
    EXPECT_EQ(matrix.getValue(1,1), 0);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_scalar_multiply, multiplicate_matrix_without_values) {
    LILMatrix matrix({{0,0},{0,0}});
    matrix.multiplicationByScalar(2);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_scalar_multiply, multiplicate_matrix_by_zero) {
    LILMatrix matrix({{1,2},{3,4}});
    matrix.multiplicationByScalar(0);
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_diagonal_multiply, multiplicate_matrix_with_full_values) {
    LILMatrix matrix({{1,2},{3,4}});
    matrix.multiplicationByDiagonalMatrix({5,6});
    EXPECT_EQ(matrix.getValue(0,0), 5);
    EXPECT_EQ(matrix.getValue(0,1), 12);
    EXPECT_EQ(matrix.getValue(1,0), 15);
    EXPECT_EQ(matrix.getValue(1,1), 24);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 4);
}
TEST(Matrix_diagonal_multiply, multiplicate_matrix_with_some_values) {
    LILMatrix matrix({{1,0},{0,4}});
    matrix.multiplicationByDiagonalMatrix({5,6});
    EXPECT_EQ(matrix.getValue(0,0), 5);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 24);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_diagonal_multiply, multiplicate_matrix_without_values) {
    LILMatrix matrix({{0,0},{0,0}});
    matrix.multiplicationByDiagonalMatrix({5,6});
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_diagonal_multiply, multiplicate_matrix_by_some_zeros) {
    LILMatrix matrix({{1,2},{3,4}});
    matrix.multiplicationByDiagonalMatrix({0,6});
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 12);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 24);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_diagonal_multiply, multiplicate_matrix_by_null_diagonal) {
    LILMatrix matrix({{1,2},{3,4}});
    matrix.multiplicationByDiagonalMatrix({0,0});
    EXPECT_EQ(matrix.getValue(0,0), 0);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 0);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 0);
}
TEST(Matrix_identity_substract, substract_full_matrix) {
    LILMatrix matrix({{0,2},{0,0}});
    matrix.identitySubtractSelf();
    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), -2);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 3);
}
TEST(Matrix_identity_substract, substract_matrix_with_some_zeros) {
    LILMatrix matrix({{0,1},{4,0}});
    matrix.identitySubtractSelf();
    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), -1);
    EXPECT_EQ(matrix.getValue(1,0), -4);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 4);
}
TEST(Matrix_identity_substract, substract_null_matrix) {
    LILMatrix matrix({{0,0},{0,0}});
    matrix.identitySubtractSelf();
    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);
}
TEST(Matrix_gaus, gaus_with_already_triangulated_matrix) {
    LILMatrix matrix({
        {1,0},
        {0,1}
    });
    vector<double> b = {5, 9};

    matrix.gaussianElimination(b);

    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);

    EXPECT_EQ(b[0], 5);
    EXPECT_EQ(b[1], 9);
}
TEST(Matrix_gaus, gaus_simple_triangulation) {
    LILMatrix matrix({
        {1,0},
        {1,1}
    });
    vector<double> b = {5, 9};

    matrix.gaussianElimination(b);

    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getRows(), 2);
    EXPECT_EQ(matrix.getCols(), 2);
    EXPECT_EQ(matrix.nzElems(), 2);

    EXPECT_EQ(b[0], 5);
    EXPECT_EQ(b[1], 4);
}
TEST(Matrix_gaus, gaus_complete_triangulation) {
    LILMatrix matrix({
        {1,2,3},
        {3,2,1},
        {2,3,1}
    });
    vector<double> b = {5, 9, 14.5};

    matrix.gaussianElimination(b);

    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 2);
    EXPECT_EQ(matrix.getValue(0,2), 3);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), -4);
    EXPECT_EQ(matrix.getValue(1,2), -8);
    EXPECT_EQ(matrix.getValue(2,0), 0);
    EXPECT_EQ(matrix.getValue(2,1), 0);
    EXPECT_EQ(matrix.getValue(2,2), -3);
    EXPECT_EQ(matrix.getRows(), 3);
    EXPECT_EQ(matrix.getCols(), 3);
    EXPECT_EQ(matrix.nzElems(), 6);

    EXPECT_EQ(b[0], 0);
    EXPECT_EQ(b[1], 5.5);
    EXPECT_EQ(b[2], -2);
}
TEST(Matrix_gaus, gaus_with_matrix_with_linearly_dependent_rows) {
    LILMatrix matrix({
        {1,2,3},
        {3,2,1},
        {1,2,3}
    });
    vector<double> b = {5, 9, 11};

    matrix.gaussianElimination(b);

    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 2);
    EXPECT_EQ(matrix.getValue(0,2), 3);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), -4);
    EXPECT_EQ(matrix.getValue(1,2), -8);
    EXPECT_EQ(matrix.getValue(2,0), 0);
    EXPECT_EQ(matrix.getValue(2,1), 0);
    EXPECT_EQ(matrix.getValue(2,2), 0);
    EXPECT_EQ(matrix.getRows(), 3);
    EXPECT_EQ(matrix.getCols(), 3);
    EXPECT_EQ(matrix.nzElems(), 5);

/*
    // NO EXISTE SOLUCIÓN
    EXPECT_EQ(b[0], 5);
    EXPECT_EQ(b[1], -6);
    EXPECT_EQ(b[2], 19/2);
    */
}

TEST(Matrix_gaus, gaus_with_4x4_matrix) {
    LILMatrix matrix({
         {1,2,3,4},
         {4,3,2,1},
         {4,3,1,2},
         {3,2,1,4}
    });

    vector<double> b = {1,1,1,1};

    matrix.gaussianElimination(b);
    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 2);
    EXPECT_EQ(matrix.getValue(0,2), 3);
    EXPECT_EQ(matrix.getValue(0,3), 4);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), -5);
    EXPECT_EQ(matrix.getValue(1,2), -10);
    EXPECT_EQ(matrix.getValue(1,3), -15);
    EXPECT_EQ(matrix.getValue(2,0), 0);
    EXPECT_EQ(matrix.getValue(2,1), 0);
    EXPECT_EQ(matrix.getValue(2,2), -1);
    EXPECT_EQ(matrix.getValue(2,3), 1);
    EXPECT_EQ(matrix.getValue(3,0), 0);
    EXPECT_EQ(matrix.getValue(3,1), 0);
    EXPECT_EQ(matrix.getValue(3,2), 0);
    EXPECT_EQ(matrix.getValue(3,3), 4);

    EXPECT_EQ(matrix.getRows(), 4);
    EXPECT_EQ(matrix.getCols(), 4);
    EXPECT_EQ(matrix.nzElems(), 10);
}

TEST(Matrix_gaus, gaus_with_4x4_matrix_v2) {
    LILMatrix matrix({
         {1,1,0,0},
         {0,1,0,0},
         {0,1,1,0},
         {2,0,1,1}
    });

    vector<double> b = {1,1,1,1};

    matrix.gaussianElimination(b);
    EXPECT_EQ(matrix.getValue(0,0), 1);
    EXPECT_EQ(matrix.getValue(0,1), 1);
    EXPECT_EQ(matrix.getValue(0,2), 0);
    EXPECT_EQ(matrix.getValue(0,3), 0);
    EXPECT_EQ(matrix.getValue(1,0), 0);
    EXPECT_EQ(matrix.getValue(1,1), 1);
    EXPECT_EQ(matrix.getValue(1,2), 0);
    EXPECT_EQ(matrix.getValue(1,3), 0);
    EXPECT_EQ(matrix.getValue(2,0), 0);
    EXPECT_EQ(matrix.getValue(2,1), 0);
    EXPECT_EQ(matrix.getValue(2,2), 1);
    EXPECT_EQ(matrix.getValue(2,3), 0);
    EXPECT_EQ(matrix.getValue(3,0), 0);
    EXPECT_EQ(matrix.getValue(3,1), 0);
    EXPECT_EQ(matrix.getValue(3,2), 0);
    EXPECT_EQ(matrix.getValue(3,3), 1);

    EXPECT_EQ(matrix.getRows(), 4);
    EXPECT_EQ(matrix.getCols(), 4);
    EXPECT_EQ(matrix.nzElems(), 5);
}