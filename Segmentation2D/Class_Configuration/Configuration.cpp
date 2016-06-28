/* ------------- Class Configuration ---------------
	file Configuration.cpp  Modified on 09/07/2012
	source file
---------------------------------------------*/

#include "../Class_Ellips/Ellips.h"
#include "Configuration.h"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "../other_functions/other_functions.h"
#include <string>
#include <omp.h>

using namespace std;

Configuration::Configuration()
{
	config = new Ellips[1];
	position = new int[1];
	data_fit = new double[1];

	config[0]= Ellips();
	position[0]=1;
	data_fit[0]=1;
	nb_Ellipses=1;
	size=1;
}

Configuration::Configuration(const Ellips & ell,int pos, double fit, int size_tot)
{
	config = new Ellips[size_tot];
	position = new int[size_tot];
	data_fit = new double[size_tot];

	config[0]=ell;
	position[0]=pos;
	data_fit[0]=fit;
	nb_Ellipses=1;
	size=size_tot;

}

Configuration::Configuration(double a_min, double a_max, int size_x, int size_y,int nb_ell, int nb_dont_accepted, double * img, double d)
{
	config = new Ellips[nb_ell];
	position = new int[nb_ell];
	data_fit = new double[nb_ell];

	int dont_accepted=0; // number of Ellipses which aren't accepted
	int inc=0; // number of Ellipses accepted

	bool inter, t_inter; // result of the intersection of 2 Ellipses
	int ind,pos, t_ind, t_inc;

	// grid 
	int nb_rows = ceil(size_y/a_max);
	int nb_col = ceil(size_x/a_max);

	Ellips new_ell;

	while ((inc<nb_ell) & (dont_accepted<nb_dont_accepted)){
		// generation of a new Ellipse
		Ellips new_ell(a_min,a_max,size_x,size_y);
		pos=min(nb_rows-1,floor(new_ell.get_cy()/a_max))*nb_col+max(1,ceil(new_ell.get_cx()/a_max));

		inter=false;
		ind=0;
		while ((inter==false) & (ind < inc)){
			if (is_neighbor(pos,position[ind],nb_rows,nb_col,a_max)){
				inter=intersect(config[ind],new_ell);
			}
			ind++;
		}

		if (!(inter)){ // we keep the Ellipse
			{
				config[inc]=new_ell;
				position[inc]=pos;
				inc++;
				dont_accepted=0;
			}
		}
		else {
			{
				dont_accepted++;
			}
		} // if (!(inter))
	} // while ((inc<50) & (dont_accepted<10))

	#pragma omp parallel for
	for(int i = 0; i < inc; i++){
		data_fit[i] = config[i].data_fiting(img,size_x,size_y,d);
	}

	nb_Ellipses=inc;
	size=inc;
}

Configuration::Configuration(double a_min, double a_max, int size_x, int size_y,int nb_ell, int nb_dont_accepted, double * grad_x, double * grad_y, double step, double epsilon, double d)
{
	config = new Ellips[nb_ell];
	position = new int[nb_ell];
	data_fit = new double[nb_ell];

	int dont_accepted=0; // number of Ellipses which aren't accepted
	int inc=0; // number of Ellipses accepted

	bool inter; // result of the intersection of 2 Ellipses
	int ind,pos;

	// grid 
	int nb_rows = ceil(size_y/a_max);
	int nb_col = ceil(size_x/a_max);

	while ((inc<nb_ell) & (dont_accepted<nb_dont_accepted)){
		// generation of a new Ellipse
		Ellips new_ell(a_min,a_max,size_x,size_y);
		pos=min(nb_rows-1,floor(new_ell.get_cy()/a_max))*nb_col+max(1,ceil(new_ell.get_cx()/a_max));

		inter=false;
		ind=0;
		while ((inter==false) & (ind < inc)){
			if (is_neighbor(pos,position[ind],nb_rows,nb_col,a_max)){
				inter=intersect(config[ind],new_ell);
			}
			ind++;
		}
		if (!(inter)){ // we keep the Ellipse
			config[inc]=new_ell;
			position[inc]=pos;

			inc++;
			dont_accepted=0;
		}
		else {
			dont_accepted++;
		} // if (ind==inc)
	} // while ((inc<50) & (dont_accepted<10))

	#pragma omp parallel for
	for(int i = 0; i < inc; i++){
		data_fit[i] = config[i].data_fiting(step,grad_x,grad_y,epsilon,size_x,size_y,d);
	}

	nb_Ellipses=inc;
	size=nb_ell;
}

Configuration::Configuration( Configuration const & other)
{
	nb_Ellipses=other.nb_Ellipses;
	size=other.size;
	config = new Ellips[size];
	position = new int[size];
	data_fit = new double[size];

	for (int i=0 ; i<nb_Ellipses ; i++){
		config[i]=other.config[i];
		position[i]=other.position[i];
		data_fit[i]=other.data_fit[i];
	}

}

Configuration::~Configuration()
{
	delete [] config;
	delete [] data_fit;
	delete [] position;
}

Configuration & Configuration::operator = ( const Configuration & other)
{
	nb_Ellipses=other.nb_Ellipses;
	size=other.size;
	delete [] config;
	delete [] data_fit;
	delete [] position;
	config = new Ellips[size];
	position = new int[size];
	data_fit = new double[size];

	for (int i=0 ; i<nb_Ellipses ; i++){
		config[i]=other.config[i];
		position[i]=other.position[i];
		data_fit[i]=other.data_fit[i];
	}

	return(*this);
}

Ellips & Configuration::get_Ellips(int i)
{
	return config[i];
}

int Configuration::get_position(int i)
{
	return position[i];
}

double Configuration::get_data_fit(int i)
{
	return data_fit[i];
}

int Configuration::get_nb_Ellipses()
{
	return nb_Ellipses;
}

void Configuration::add_Ellips(const Ellips & ell,int pos, double fit)
{
	nb_Ellipses++;
	if (nb_Ellipses>size){
		delete [] config;
		delete [] data_fit;
		delete [] position;
	
		Configuration temp(*this);
		config = new Ellips[nb_Ellipses];
		position = new int[nb_Ellipses];
		data_fit = new double[nb_Ellipses];

		for (int i=0 ; i<nb_Ellipses-1 ; i++){
			config[i]=temp.config[i];
			position[i]=temp.position[i];
			data_fit[i]=temp.data_fit[i];
		}
	}

	config[nb_Ellipses-1]=ell;
	position[nb_Ellipses-1]=pos;
	data_fit[nb_Ellipses-1]=fit;

}

void Configuration::save_config(string file_name)
{
	ofstream file(file_name.c_str(), ios::out | ios::trunc); // To open the file
	if(file) {
		file << nb_Ellipses << endl;
		for (int i = 0 ; i < nb_Ellipses ; i++){
			file << config[i].get_a() << "\t" << config[i].get_b() << "\t" << config[i].get_theta() << "\t" << config[i].get_cx() << "\t" << config[i].get_cy() << endl;
		}
		file.close(); // To close the file
	}
	else{
		cerr << "Error to open the file " << file_name << endl;
	} // if(file)
}

