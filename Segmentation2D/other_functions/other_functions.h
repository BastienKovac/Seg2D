/*------------ Other Functions --------------
	file other_functions.h  Modified on 08/23/2012
	header file

Description :
	Contains some function which cannot be in the classes like 
	the functions min and max...
---------------------------------------------*/

#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H

#include <highgui.h>

//-- Definition of a function min
#define min(a,b) (a<=b?a:b)
//-- Definition of a function max
#define max(a,b) (a>=b?a:b)

//-- To know the number of the bax (i,j) in the grid (reshape grid => vector)
int ij_to_k (int i , int  j, double size_x, double size_y, double a_max);

//-- Inverse of the previous function (vector => grid)
void k_to_ij (int k,int & i, int & j,int nb_rows , int nb_col, double a_max);

//-- To know (in vector representation) if two box are 25 neighbor (5x5 box)
bool is_neighbor (int k1,int k2,int nb_rows , int nb_col, double a_max);

//-- Compute 1-exp(-alpha/(1-x))*exp(alpha) ([0,1] --> [0,1])
double shift_cost_exp( double x, double alpha);

//-- Compute exp(0.6931*(x-1)/(1-s)) [-1,1] --> [0,1] for all x < s the function is < to 0.5 (s is an acceptance limit for the graph)
double shift_cost_exp1( double x, double s);

//-- Component x of the gradient of an image
//	    img : array which contains the image
//	    grad : array to put the result (same size as img)
//		size_x size_y : dimensions of the image
void grad_X (double * img, double * grad, int size_x, int size_y);

//-- Component y of the gradient of an image
//	    img : array which contains the image
//	    grad : array to put the result (same size as img)
//		size_x size_y : dimensions of the image
void grad_Y (double * img, double * grad, int size_x, int size_y);

//-- Interpol the image
//		x,y : coordinates of the point which we want to know the value
//	    img : array which contains the image
//		size_x size_y : dimensions of the image
double interp_grad(double x, double y, double * img, int size_x);

//-- To know the min and max value of an image img
//		img : image of type IplImage
//      min_val max_val : variable to put the result
void min_max_val (IplImage* img, double & min_val, double & max_val);

//-- To know the min and max value of an array mat
//		mat : array which contains the image
//      min_val max_val : variable to put the result
//		size_x size_y : dimensions of the image
void min_max_val (double * mat, double & min_val, double & max_val, int size_x, int size_y);

//-- To convert an IplImage (array of uchar) into an array of double
void convert_char_to_double(IplImage* img, double * mat);

//-- To convert an array of double into  an IplImage (array of uchar)
//		need to know the max value of the IplImage codage (ex 255 for image 8bit...)
void convert_double_to_char(IplImage* img, double * mat,double max_val);

//-- To smooth an image
//		convolution with a gaussian kernel (variance sigma)
void smooth_img(double * img, double sigma, int size_x, int size_y);

#endif
