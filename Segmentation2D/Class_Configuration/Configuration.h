/*------------- Class Configuration ---------------
	file Configuration.h  Modified on 09/07/2012
	header file

Description :
	Class of Configuration of Ellipses define by 3 arrays :
		- config : array to save the Ellipses of the configuration
		- position : array to save the position (in the partition of the image with squares)
		- data_fit : array to save the the data fitting of the Ellipses
---------------------------------------------*/

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "../Class_Ellips/Ellips.h"
#include <string>

using namespace std;

class Configuration
{
public :

	//////////////////////////////////////////
	///////////// CONSTRUCTORS ///////////////
	//////////////////////////////////////////

	//-- Default constructor (configuration of one Ellipse with position 1 and data fit 1)
	Configuration();

	//-- Second constructor (with parameters)
	// Generate a configuration of one Ellipse given as parameter but defines the array with a size equal to
	// the parameter size_tot.
	Configuration(const Ellips & ell,int pos, double fit, int size_tot);

	//-- Third constructor (with parameters)
	// Generate a random configuration of Ellipses without intersections
	//		a_min : minimum value accepted for the axis
	//		a_max : maximum value accepted for the axis
	//		size_x, size_y : size of the image where the Ellipses live
	//		nb_ell : number of Ellipses wanted (size of the arrays config, position and data_fit)
	//      nb_dont_accepted : number of Ellipse don't accepted before to stop the algorithm
	//						   ( example : if we want 50 Ellipses (nb_ell = 50) and because of the intersection constraint,
	//							the algorithm rejects nb_dont_accepted Ellipses successively the algorithm is stopped and return
	//							the real number of Ellipses generated )
	//		img : image in which we want to find objects
	//		d : acceptance threshold
	//-- Version 1 -- (with the function data fit 1)
	Configuration(double a_min, double a_max, int size_x, int size_y,int nb_ell, int nb_dont_accepted,double * img, double d);
	//-- Version 2 -- (with the function data fit 2)
	//		grad_x : component x of the gradient of the image
	//		grad_y : component y of the gradient of the image
	//		epsilon : regularisation constant
	//		step : step for the numerical integration
	Configuration(double a_min, double a_max, int size_x, int size_y,int nb_ell, int nb_dont_accepted,double* grad_x, double* grad_y, double step, double epsilon, double d);

	//-- Copy constructor
	Configuration( Configuration const & other);

	//-- Destructor
	~Configuration();

	//////////////////////////////////////////
	/////////////// OPERATORS ////////////////
	//////////////////////////////////////////

	//-- Overload of the = operator
	Configuration & operator = ( const Configuration & other);

	///////////////////////////////////////////
	/////////////// ACCESSORS /////////////////
	///////////////////////////////////////////

	//-- To get the Ellipse at the index i in the array config
	Ellips & get_Ellips(int i);

	//-- To get the position at the index i in the array position
	int get_position(int i);

	//-- To get the data_fit at the index i in the array data_fit
	double get_data_fit(int i);

	//-- To get the number of Ellipses in the configuration
	int get_nb_Ellipses();

	///////////////////////////////////////////
	///////////////// METHODS /////////////////
	///////////////////////////////////////////

	//-- To add an Ellipse to the configuration
	void add_Ellips(const Ellips & ell,int pos, double fit);

	//-- To save the parameters of the Ellipses into a file
	void save_config(string file_name);

private :

	int nb_Ellipses, size; // number of Ellipses and size of the arrays (size must be greater than nb_Ellipss);
	Ellips * config ;
	int * position ;
	double * data_fit ;

};

#endif
