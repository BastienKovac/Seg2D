/*
 * matrix.h
 *
 *  Created on: Jul 1, 2016
 *      Author: bastienkovac
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdexcept>
#include <cmath>

class Matrix {

private:

	double* internal;
	int rows, cols;

public:

	Matrix();

	Matrix(int r, int c, double val[]);

	Matrix(int r, int c);

	Matrix( const Matrix& other);

	~Matrix();

	void set(int i, int j, double value);

	double get(int i, int j);

	Matrix operator * ( const Matrix & other);

	Matrix operator * ( const double & scalar);

	Matrix operator / ( const double & scalar);

	Matrix operator + ( const Matrix & other);

	Matrix operator - ( const Matrix & other);

	Matrix transpose();

	double norm2();

};

#endif /* MATRIX_H_ */
