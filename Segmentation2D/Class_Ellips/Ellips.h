/*------------- Class Ellips ---------------
	file Ellips.h  Modified on 08/23/2012
	header file

Description :
	Class of Ellipses objects which are defined by five parameters :
		- Semi major axis
		- Semi minor axis
		- Angle with horizontal
		- Coordinates of the center
---------------------------------------------*/

#ifndef Ellips_H
#define Ellips_H

#include "opencv/highgui.h"

class Ellips
{
public:

	//////////////////////////////////////////
	///////////// CONSTRUCTORS ///////////////
	//////////////////////////////////////////

	//-- Default constructor (create a circle of radius 1, centered at (0,0))
	Ellips();

	//-- Constructor which randomly generate an Ellipse
	//		a_min : minimum value accepted for the axis
	//		a_max : maximum value accepted for the axis
	//		size_x, size_y : size of the image (interval where the Ellipse can exist)
	Ellips(double a_min, double a_max, int size_x, int size_y);

	//-- Constructor which takes parameters of the Ellipse
	//		major_a : semi major axis
	//		minor_a : semi minor axis
	//		angle : angle with horizontal
	//		cx, cy : coordinates of the center 
	Ellips(double major_a, double minor_a, double angle, double cx, double cy);

	//-- Copy constructor
	Ellips( Ellips const & ell);

	//-- Destructor
	~Ellips();

	//////////////////////////////////////////
	/////////////// OPERATORS ////////////////
	//////////////////////////////////////////

	//-- Overload of the = operator
	Ellips & operator = ( const Ellips & other);

	///////////////////////////////////////////
	/////////////// ACCESSORS /////////////////
	///////////////////////////////////////////

	//-- To get the semi major axis
	double get_a() const;

	//-- To get the semi minor axis
	double get_b() const;

	//-- To get the coordinate x of the center
	double get_cx() const;

	//-- To get the coordinate y of the center
	double get_cy() const;

	//-- To get the angle theta
	double get_theta() const;

	///////////////////////////////////////////
	///////////////// METHODS /////////////////
	///////////////////////////////////////////

	//-- To known if (x-c)'*A*(x-c) <= level
	//		x, y : coordinates of the point which we want to know if it is inside the Ellipse
	//		level : level set which define the Ellipse
	bool inside(double x, double y, double level) const;

	//-- Compute the value x'*A*x
	//		x, y : coordinates of the point which we want to know the level set
	double level_set(double x, double y) const;

	//-- To know if 2 Ellipses intersect (friend method)
    friend bool intersect( Ellips& one,  Ellips& two);

	//-- Display the parameters of the Ellipse
	void display();

	//-- To know how an Ellipse is good or not with respect to an image
	//		im : Image
	//		this function return a number between 0 and 1
	//			- near 0 --> bad
	//			- near 1 --> good
	//		d : acceptance threshold
	//-- Version 1 --
	double data_fiting (double* im ,int size_x, int size_y, double d) const;
	//-- Version 2 --
	//		step : step for the numerical integration
	//		epsilon : regularization constant
	double data_fiting (double step, double * grad_x , double* grad_y, double epsilon, int size_x, int size_y, double d) const;

private:

	double a,b,theta; // a : Semi major axis, b : Semi minor axis, theta : angle with horizontal
	double A[2][2]; // matrix of the bilinear form which represent the Ellipse
	double A_demi[2][2]; // root of A such as A_demi*A_demi=A
	double A_inv_demi[2][2]; // inverse of A_demi
	double c[2]; // vector which contains the coordinates of the center of the Ellipse

};


//-- To know if 2 Ellipses intersect
bool intersect(  Ellips & one,  Ellips & two);

#endif
