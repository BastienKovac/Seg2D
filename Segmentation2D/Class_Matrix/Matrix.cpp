/*
 * Matrix.cpp
 *
 *  Created on: Jul 1, 2016
 *      Author: bastienkovac
 */

#include "Matrix.h"

Matrix::Matrix() {

	rows = 0;
	cols = 0;

	internal = new double[rows * cols];

}

Matrix::Matrix(int r, int c, double val[]) {

	rows = r;
	cols = c;

	internal = new double[rows * cols];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			internal[i + j * cols] = val[i + j * cols];
		}
	}

}

Matrix::Matrix(int r, int c) {

	rows = r;
	cols = c;

	internal = new double[rows * cols];

}

Matrix::Matrix(const Matrix& other) {

	rows = other.rows;
	cols = other.cols;

	internal = new double[rows * cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			internal[i + j * cols] = other.internal[i + j * cols];
		}
	}

}

Matrix::~Matrix() {
	delete[] internal;
}

void Matrix::set(int i, int j, double value) {

	internal[i + j * cols] = value;

}

double Matrix::get(int i, int j) {

	return internal[i + j * cols];

}

Matrix Matrix::operator *(const Matrix & other) {

	if (!(cols == rows)) {
		throw std::runtime_error("Dimensions are wrong !");
	}

	int ret_rows = rows, ret_cols = other.cols;

	double ret_array[ret_rows * ret_cols];
	double tmp;

	for (int i = 0; i < ret_rows; i++) {
		for (int j = 0; j < ret_cols; j++) {
			tmp = 0;
			for (int k = 0; k < other.rows; k++) {
				tmp += (internal[i + k * cols]
						* other.internal[k + j * other.cols]);
			}
			ret_array[i + j * cols] = tmp;
		}
	}

	return Matrix(ret_rows, ret_cols, ret_array);

}

Matrix Matrix::operator *(const double & scalar) {

	double ret_array[rows * cols];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret_array[i + j * cols] = (internal[i + j * cols] * scalar);
		}
	}

	return Matrix(rows, cols, ret_array);

}

Matrix Matrix::operator /(const double & scalar) {

	double ret_array[rows * cols];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret_array[i + j * cols] = (internal[i + j * cols] / scalar);
		}
	}

	return Matrix(rows, cols, ret_array);

}

Matrix Matrix::operator +(const Matrix & other) {

	if (rows != other.rows || cols != other.cols) {
		throw std::runtime_error("Matrix are not the same dimension");
	}

	double ret_array[rows * cols];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret_array[i + j * cols] = (internal[i + j * cols]
					+ other.internal[i + j * other.cols]);
		}
	}

	return Matrix(rows, cols, ret_array);

}

Matrix Matrix::operator -(const Matrix & other) {

	if (rows != other.rows || cols != other.cols) {
		throw std::runtime_error("Matrix are not the same dimension");
	}

	double ret_array[rows * cols];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret_array[i + j * cols] = (internal[i + j * cols]
					- other.internal[i + j * other.cols]);
		}
	}

	return Matrix(rows, cols, ret_array);

}

Matrix Matrix::transpose() {

	double ret_array[cols * rows];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret_array[i + j * rows] = internal[j + i * cols];
		}
	}

	return Matrix(cols, rows, ret_array);

}

double Matrix::norm2() {

	double ret = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret += (internal[i + j * cols] * internal[i + j * cols]);
		}
	}

	return sqrt(ret);

}

