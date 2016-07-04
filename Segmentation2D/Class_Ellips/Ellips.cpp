/* ------------- Class Ellips ---------------
 file Ellips.cpp  Modified on 08/23/2012
 source file
 ---------------------------------------------*/

#include "Ellips.h"

/***** Note ******
 In this class we use the blas fonctions to compute matrix product.
 Representation of the matrix :
 M=| 2 3 |
 | 4 9 |   =>  M[2][2]={ {2 , 4} , { 3 , 9} } ( Rows major !!!)

 ***********/

Ellips::Ellips() {
	a = 1;
	b = 1;
	theta = 0;
	c[0] = 0;
	c[1] = 0;

	double P[4] = { cos(-theta), sin(-theta), -sin(-theta), cos(-theta) };

	//double trans_P[2][2]={{cos(-theta),-sin(-theta)},{sin(-theta),cos(-theta)}}; // transpose of P

	double D[4] = { 1 / (a * a), 0, 0, 1 / (b * b) };
	double D_demi[4] = { 1 / a, 0, 0, 1 / b }; // root of D
	double D_inv_demi[4] = { a, 0, 0, b }; // inverse of D_demi

	Matrix P_mat = Matrix(2, 2, P);
	Matrix D_mat = Matrix(2, 2, D);
	Matrix D_demi_mat = Matrix(2, 2, D_demi);
	Matrix D_inv_demi_mat = Matrix(2, 2, D_inv_demi);

	// Compute the matrix A P'*D*P
	A = P_mat.transpose() * D_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A[0][0],&n);
	 */

	// Compute the matrix A_demi P'*D_demi*P
	A_demi = P_mat.transpose() * D_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_demi[0][0],&n);
	 */

	// Compute the matrix A_inv_demi P'*D_inv_demi*P
	A_inv_demi = P_mat.transpose() * D_inv_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_inv_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_inv_demi[0][0],&n);
	 */
}

Ellips::Ellips(double major_a, double minor_a, double angle, double cx,
		double cy) {
	a = Max(major_a, minor_a);
	b = Min(major_a, minor_a);
	theta = fmod(angle, double(M_PI)); // to have a number between 0 and pi
	c[0] = cx;
	c[1] = cy;

	double P[4] = { cos(-theta), sin(-theta), -sin(-theta), cos(-theta) };

	//double trans_P[2][2]={{cos(-theta),-sin(-theta)},{sin(-theta),cos(-theta)}}; // transpose of P

	double D[4] = { 1 / (a * a), 0, 0, 1 / (b * b) };
	double D_demi[4] = { 1 / a, 0, 0, 1 / b }; // root of D
	double D_inv_demi[4] = { a, 0, 0, b }; // inverse of D_demi

	Matrix P_mat = Matrix(2, 2, P);
	Matrix D_mat = Matrix(2, 2, D);
	Matrix D_demi_mat = Matrix(2, 2, D_demi);
	Matrix D_inv_demi_mat = Matrix(2, 2, D_inv_demi);

	// Compute the matrix A P'*D*P
	A = P_mat.transpose() * D_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A[0][0],&n);
	 */

	// Compute the matrix A_demi P'*D_demi*P
	A_demi = P_mat.transpose() * D_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_demi[0][0],&n);
	 */

	// Compute the matrix A_inv_demi P'*D_inv_demi*P
	A_inv_demi = P_mat.transpose() * D_inv_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_inv_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_inv_demi[0][0],&n);
	 */
}

Ellips::Ellips(double a_min, double a_max, int size_x, int size_y) {
	double aux;
	if (a_min > a_max) {
		aux = a_min;
		a_min = a_max;
		a_max = aux;
	}

	a = (rand() / double(RAND_MAX)) * (a_max - a_min) + a_min;
	theta = (rand() / double(RAND_MAX)) * M_PI;
	c[0] = (rand() / double(RAND_MAX)) * (size_x - 1);
	b = (rand() / double(RAND_MAX)) * (a - a_min) + a_min;
	c[1] = (rand() / double(RAND_MAX)) * (size_y - 1);

	double P[4] = { cos(-theta), sin(-theta), -sin(-theta), cos(-theta) };

	//double trans_P[2][2]={{cos(-theta),-sin(-theta)},{sin(-theta),cos(-theta)}}; // transpose of P

	double D[4] = { 1 / (a * a), 0, 0, 1 / (b * b) };
	double D_demi[4] = { 1 / a, 0, 0, 1 / b };
	; // root of D
	double D_inv_demi[4] = { a, 0, 0, b }; // inverse of D_demi

	Matrix P_mat = Matrix(2, 2, P);
	Matrix D_mat = Matrix(2, 2, D);
	Matrix D_demi_mat = Matrix(2, 2, D_demi);
	Matrix D_inv_demi_mat = Matrix(2, 2, D_inv_demi);

	// Compute the matrix A P'*D*P
	A = P_mat.transpose() * D_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A[0][0],&n);
	 */

	// Compute the matrix A_demi P'*D_demi*P
	A_demi = P_mat.transpose() * D_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_demi[0][0],&n);
	 */

	// Compute the matrix A_inv_demi P'*D_inv_demi*P
	A_inv_demi = P_mat.transpose() * D_inv_demi_mat * P_mat;

	/*
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&D_inv_demi[0][0],&n,&P[0][0],&n,&beta,&temp[0][0],&n);
	 F77_NAME(dgemm)(&trans,&trans,&n,&n,&n,&alpha,&trans_P[0][0],&n,&temp[0][0],&n,&beta,&A_inv_demi[0][0],&n);
	 */
}

Ellips::Ellips(Ellips const & ell) {
	a = ell.a;
	b = ell.b;
	theta = ell.theta;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			/*
			 A[i][j]=ell.A[i][j];
			 A_demi[i][j]=ell.A_demi[i][j];
			 A_inv_demi[i][j]=ell.A_inv_demi[i][j];
			 */
			A = Matrix(ell.A);
			A_demi = Matrix(ell.A_demi);
			A_inv_demi = Matrix(ell.A_inv_demi);
		}
		c[i] = ell.c[i];
	}
}

Ellips::~Ellips() {
}

Ellips & Ellips::operator =(const Ellips & other) {
	a = other.a;
	b = other.b;
	theta = other.theta;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			A = Matrix(other.A);
			A_demi = Matrix(other.A_demi);
			A_inv_demi = Matrix(other.A_inv_demi);
		}
		c[i] = other.c[i];
	}

	return (*this);
}

double Ellips::get_a() const {
	return a;
}

double Ellips::get_b() const {
	return b;
}

double Ellips::get_cx() const {
	return c[0];
}

double Ellips::get_cy() const {
	return c[1];
}

double Ellips::get_theta() const {
	return theta;
}

bool Ellips::inside(double x, double y, double level) const {
	double result;
	//double temp1[2]={x-c[0],y-c[1]}; // temp1 <-- ([x;y] - c)

	Matrix tmp = Matrix(2, 1);
	tmp.set(0, 0, x - c[0]);
	tmp.set(1, 0, y - c[1]);

	result = (tmp.transpose() * A * tmp).get(0, 0);

	/*
	 F77_NAME(dgemv)(&trans,&n,&n,&alpha,&A[0][0],&n,temp1,&inc,&beta,temp2,&inc); // temp2 <-- A*temp1 = A*([x;y] - c)
	 result=F77_NAME(ddot)(&n, temp1, &inc,temp2, &inc); //result <-- temp1.temp2 = ([x;y] - c)'*A*([x;y] - c)
	 */

	return (result <= level);

}

double Ellips::level_set(double x, double y) const {
	double result;
	//double temp1[2]={x,y}; // temp1 <-- [x;y]

	Matrix tmp = Matrix(2, 1);
	tmp.set(0, 0, x);
	tmp.set(1, 0, y);

	result = (tmp.transpose() * A * tmp).get(0, 0);

	/*
	 F77_NAME(dgemv)(&trans,&n,&n,&alpha,&A[0][0],&n,temp1,&inc,&beta,temp2,&inc); // temp2 <-- A*temp1 = A*[x;y]
	 result=F77_NAME(ddot)(&n, temp1, &inc,temp2, &inc); //result <-- temp1.temp2 = [x;y]'*A*[x;y]
	 */

	return result;
}

bool intersect(Ellips& one, Ellips& two) {
	//--- lipschitz constants -
	double L1 = 2 * two.a * two.a / (one.b * one.b); // definition of the lipschitz constants
	double mu = 2 * two.b * two.b / (one.a * one.a);

	//double x[2]={two.c[0]-one.c[0],two.c[1]-one.c[1]}; // x <-- (c2 - c1)

	Matrix tmp = Matrix(2, 1);
	tmp.set(0, 0, two.c[0] - one.c[0]);
	tmp.set(1, 0, two.c[1] - one.c[1]);

	double cond1 = (tmp.transpose() * one.A * tmp).get(0, 0);

	/*
	 F77_NAME(dgemv)(&trans,&n,&n,&alpha,&one.A[0][0],&n,x,&inc,&beta,res,&inc); // res <-- A1*x = A1*(c2 - c1)
	 double cond1= F77_NAME(ddot)(&n, x, &inc,res, &inc); // cond1 <-- x.res = (c2 - c1)'*A1*(c2 - c1)
	 */

	tmp.set(0, 0, one.c[0] - two.c[0]);
	tmp.set(1, 0, one.c[1] - two.c[1]);

	double cond2 = (tmp.transpose() * two.A * tmp).get(0, 0);

	/*
	 F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A[0][0],&n,x,&inc,&beta,res,&inc); // res <-- A2*x = A2*(c1 - c2)
	 double cond2= F77_NAME(ddot)(&n, x, &inc,res, &inc); // cond2 <-- x.res = (c1 - c2)'*A2*(c1 - c2)
	 */

	if ((cond1 <= 1) | (cond2 <= 1)) {
		return true; // the center of one Ellips is in the other Ellips => intersect
	} else {
		// Method of projected gradient
		//double y[2]={(2*one.get_cx()+two.get_cx())/3,(2*one.get_cy()+two.get_cy())/3};
		bool cond = true;
		double norm;
		double lambda_min;

		Matrix y = Matrix(2, 1);
		y.set(0, 0, (2 * one.get_cx() + two.get_cx()) / 3);
		y.set(1, 0, (2 * one.get_cy() + two.get_cy()) / 3);

		Matrix c1 = Matrix(2, 1);
		c1.set(0, 0, one.c[0]);
		c1.set(1, 0, one.c[1]);

		Matrix c2 = Matrix(2, 1);
		c2.set(0, 0, two.c[0]);
		c2.set(1, 0, two.c[1]);

		// Computation of the gradient [grad=2*A2_inv_demi*A1*(A2_inv_demi*y-c1)]

		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_inv_demi[0][0],&n,y,&inc,&beta,res,&inc); //res <-- A2_inv_demi*y
		Matrix res = two.A_inv_demi * y;

		//F77_NAME(daxpy)(&n, &alpha ,one.c, &inc,res, &inc); // res <-- res - c1 = A2_inv_demi*y -c1
		res = res - c1;

		//F77_NAME(dscal)(&n, &alpha, res,&inc); // res <-- 2*res = 2*(A2_inv_demi*y -c1)
		res = res * 2;

		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&one.A[0][0],&n,res,&inc,&beta,grad,&inc); //grad <-- A1*res= 2*A1*(A2_inv_demi*y -c1)
		//F77_NAME(dcopy)(&n, grad,&inc,res,&inc); // res <-- grad
		Matrix grad = one.A * res;

		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_inv_demi[0][0],&n,res,&inc,&beta,grad,&inc); // grad <-- A2_inv_demi*res = 2*A2_inv_demi*A1*(A2_inv_demi*y -c1)
		grad = two.A_inv_demi * res;

		while (cond) {

			//F77_NAME(dcopy)(&n, y,&inc,y_old,&inc); // y_old <-- y
			Matrix y_old = Matrix(y);

			// Gradient descent

			//F77_NAME(daxpy)(&n, &alpha,grad, &inc,y, &inc); // y=y-grad*2/(L1+mu)
			y = y - (grad * (2 / (L1 + mu)));

			//F77_NAME(dcopy)(&n, y,&inc,res,&inc); // res <-- y
			res = Matrix(y);

			// Projection of y on the Ellips two
			//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_demi[0][0],&n,two.c,&inc,&beta,res,&inc); //res <--  res - A2_demi*c2 = y -A2_demi*c2
			//norm=F77_NAME(dnrm2)(&n,res, &inc); // norm(res) = norm(y -A2_demi*c2)
			res = res - (two.A_demi * c2);
			norm = res.norm2();

			if (norm > 1) {
				//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_demi[0][0],&n,two.c,&inc,&beta,res,&inc); // y <-- A2_demi*c2 + res/norm = A2_demi*c2 + (y -A2_demi*c2)/norm(y -A2_demi*c2)
				//F77_NAME(dcopy)(&n, res,&inc,y,&inc); // y <-- res
				y = (two.A_demi * c2) + (res / norm);

			}
			// Computation of the gradient
			//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_inv_demi[0][0],&n,y,&inc,&beta,res,&inc); //res <-- A2_inv_demi*y
			res = two.A_inv_demi * y;

			//F77_NAME(daxpy)(&n, &alpha ,one.c, &inc,res, &inc); // res <-- res - c1 = A2_inv_demi*y -c1
			res = res - c1;

			//F77_NAME(dscal)(&n, &alpha, res,&inc); // res <-- 2*res = 2*(A2_inv_demi*y -c1)
			res = res * 2;

			//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&one.A[0][0],&n,res,&inc,&beta,grad,&inc); //grad <-- A1*res= 2*A1*(A2_inv_demi*y -c1)
			grad = one.A * res;

			//F77_NAME(dcopy)(&n, grad,&inc,res,&inc); // res <-- grad
			res = grad;

			//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_inv_demi[0][0],&n,res,&inc,&beta,grad,&inc); // grad <-- A2_inv_demi*res = 2*A2_inv_demi*A1*(A2_inv_demi*y -c1)
			grad = two.A_inv_demi * res;

			// Update the condition cond
			//F77_NAME(dcopy)(&n, y,&inc,res,&inc); // res <-- y
			res = y;

			//F77_NAME(daxpy)(&n, &alpha ,y_old, &inc,res, &inc); // res <-- res - y_old = y - y_old
			res = res - y_old;

			//norm=F77_NAME(dnrm2)(&n,res, &inc); // norm(res) = norm(y - y_old)
			res.norm2();

			cond = (norm > 0.1);
		} // while (cond)

		// Value of the minimal level set of the Ellips one which intersect with the Ellips two
		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&two.A_inv_demi[0][0],&n,y,&inc,&beta,res,&inc); // res <-- A2_inv_demi*y
		res = two.A_inv_demi * y;

		//F77_NAME(daxpy)(&n, &alpha ,one.c, &inc,res, &inc); // res <-- res - c1
		res = res - c1;

		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&one.A[0][0],&n,res,&inc,&beta,y,&inc); // y <-- A1*res = A1*(A2_inv_demi*y - c1)
		y = one.A * res;

		//lambda_min=F77_NAME(ddot)(&n, res, &inc,y, &inc); // lambda_min = y.res = (A2_inv_demi*y - c1)'*A1*(A2_inv_demi*y - c1)
		lambda_min = (res.transpose() * y * res).get(0, 0);

		lambda_min = lambda_min - 1;

		if (lambda_min < 0) {
			return true;
		} else {
			return false;
		}
	} // if ((cond1 <= 1 )| (cond2 <= 1)) -> else
}

void Ellips::display() {
	std::cout << "Semi major axe : " << a << std::endl;
	std::cout << "Semi minor axe : " << b << std::endl;
	std::cout << "Angle : " << theta << std::endl;
	std::cout << "Coordinates of the center, x : " << c[0] << " y : " << c[1]
			<< std::endl;
}

double Ellips::data_fiting(double * im, int size_x, int size_y,
		double d) const {
	double rho = 2;
	Ellips bound(a + rho, b + rho, theta, c[0], c[1]);

	double mu_in = 0;
	double mu_out = 0;
	double q;
	int nb_in = 0;
	int nb_out = 0;
	int j_stop = Min(size_y - 1, ceil(c[1] + (a + rho)));
	int i_stop = Min(size_x - 1, ceil(c[0] + (a + rho)));

	// it first passes through the y because of the structure of the memory : the x are contiguous
#pragma omp parallel for reduction(+:mu_in,mu_out,nb_in,nb_out) collapse(2) if(using_multithread)
	for (int j = Max(0, floor(c[1] - (a + rho))); j <= j_stop; j++) {
		for (int i = Max(0, floor(c[0] - (a + rho))); i <= i_stop; i++) {
			if (bound.inside(i, j, 1)) {
				if (inside(i, j, 1)) {
					mu_in += im[i + j * size_x];
					nb_in++;
				} else {
					mu_out += im[i + j * size_x];
					nb_out++;
				}
			}
		}
	}

	mu_in = mu_in / double(nb_in);
	mu_out = mu_out / double(nb_out);

	q = mu_in - mu_out;
	return shift_cost_exp1(q, d);
}

double Ellips::data_fiting(double step, double* grad_x, double* grad_y,
		double epsilon, int size_x, int size_y, double d) const {
	// Parameter of the Ellipse
	double per = 0;

	// integration of the image on the Ellipse
	double u = 0;

	// discretisation [0:2pi]
	int nb = int(2 * M_PI / step);

	double r, th = 0;
	double x, y, norme, ix, iy;
	double nor[2], res[2];

	// Parameters for the BLAS functions
	int n = 2;

	double ctheta = cos(theta);
	double stheta = sin(theta);

#pragma omp parallel for if(using_multithread)
	for (int i = 0; i < nb; i++) {
		double cth = cos(th);
		double sth = sin(th);

		r = sqrt(a * a * cth * cth + b * b * sth * sth);
		x = c[0] + a * cth * ctheta - b * sth * stheta;
		y = c[1] + a * cth * stheta + b * sth * ctheta;
		// Compute the normal to the Ellipse at the point (x,y)
		res[0] = x - c[0];
		res[1] = y - c[1];
		//F77_NAME(dgemv)(&trans,&n,&n,&alpha,&A[0][0],&n,res,&inc,&beta,nor,&inc);
		res[0] *= n;
		res[1] *= n;

		norme = sqrt(nor[0] * nor[0] + nor[1] * nor[1]);
		nor[0] = nor[0] / norme;
		nor[1] = nor[1] / norme;
		if ((x >= 0) && (x < size_x - 1) && (y >= 0) && (y < size_y - 1)) {
			ix = interp_grad(x, y, grad_x, size_x);
			iy = interp_grad(x, y, grad_y, size_x);
			norme = sqrt(ix * ix + iy * iy);
			u -= r * step
					* ((ix * nor[0] + iy * nor[1])
							/ sqrt(norme * norme + epsilon * epsilon));
			per += r * step;
		}
		th += step;
	}

	return shift_cost_exp1(u / per, d);
}
