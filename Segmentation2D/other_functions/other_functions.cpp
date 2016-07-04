/* ------------ Other Functions --------------
	file other_functions.cpp  Modified on 08/03/2012 
	source file
---------------------------------------------*/

#include "other_functions.h"

int ij_to_k (int i , int  j, double size_x, double size_y, double a_max)
{
	int nb_rows = ceil(size_y/a_max);
	int nb_col = ceil(size_x/a_max);

	if (i> nb_rows){
		std::cout << " Index i out of the image " << std::endl;
		return 0;
	}
	else{ if (j> nb_col){
		std::cout << " Index j out of the image" << std::endl;
		return 0;
		}
		else {
			return (i-1)*nb_col+j;
		}
	}
}

void k_to_ij (int k,int & i, int & j,int nb_rows , int nb_col, double a_max)
{
	if (k > nb_rows*nb_col){
		std::cout << "Index k out of the image " << k << std::endl;
	}
	else {
		i=ceil(float(k)/float(nb_col));
		j=(fmod(float(k),float(nb_col))==0) ? nb_col : fmod(float(k),float(nb_col)) ;
	}
}

bool is_neighbor (int k1,int k2,int nb_rows , int nb_col, double a_max)
{
	int i1,j1,i2,j2;
	k_to_ij(k1,i1,j1,nb_rows,nb_col,a_max);
	k_to_ij(k2,i2,j2,nb_rows,nb_col,a_max);

	if ((i2 >= Max(1,i1-2)) & (i2 <= Min(nb_rows,i1+2)) & (j2 >= Max(1,j1-2)) & (j2 <= Min(nb_col,j1+2))){
		return true;
	}
	else{
		return false;
	}

}

double shift_cost_exp( double x, double alpha)
{
	return x==1 ?1 : 1-exp(-alpha/(1-x))*exp(alpha);
}

double shift_cost_exp1( double x, double s)
{
	return exp(0.6931*(x-1)/(1-s));
}

void grad_X (double * img, double * grad, int size_x, int size_y)
{

	for (int y = 0; y < size_y; y++)
	{
		for (int x = 0; x < size_x-1 ; x++)
		{
			grad[x+y*size_x]= img[x+y*size_x+1] -  img[x+y*size_x];
		}
		grad[size_x-1+y*size_x]=img[size_x-1+y*size_x];
	}
}

void grad_Y (double * img, double * grad, int size_x, int size_y)
{
	for (int x = 0; x < size_x; x++)
	{
		for (int y = 0; y < size_y-1; y++)
		{
			grad[x+y*size_x]= img[x+y*size_x+size_x] -  img[x+y*size_x];
		}
		grad[x+(size_y-1)*size_x]=img[x+(size_y-1)*size_x];
	}

}

double interp_grad(double x, double y, double * img, int size_x)
{
	int j = floor(y);
	int i = floor(x);
	double u;
	u=img[i+j*size_x]*(j+1-y)*(i+1-x)+img[i+j*size_x+1]*(x-i)*(j+1-y)+img[i+j*size_x+size_x]*(i+1-x)*(y-j)+img[i+j*size_x+size_x+1]*(x-i)*(y-j);
	return u;

}

void min_max_val (IplImage* img, double & min_val, double & max_val)
{
	uchar *p;
	p = cvPtr2D (img, 0, 0, NULL);
	min_val=*p; max_val=*p;
	for (int y = 0; y < img->height; y++)
	{
		for (int x = 0; x < img->width ; x++)
		{
			// pointer on the pixel of the images
			p = cvPtr2D (img, y, x, NULL);
			min_val=Min(min_val,(double)(*p));
			max_val=Max(max_val,(double)(*p));
		}
	}
}

void min_max_val (double * mat, double & min_val, double & max_val, int size_x, int size_y)
{
	max_val=mat[0];
	min_val=mat[0];
	for (int y=0 ; y < size_y ; y++){
		for (int x =0 ; x < size_x ; x++){
			min_val=Min(min_val,mat[x+y*size_x]);
			max_val=Max(max_val,mat[x+y*size_x]);
		}
	}

}

void convert_char_to_double(IplImage* img, double * mat)
{
	uchar *p;
	double min_val, max_val;
	min_max_val(img,min_val,max_val);
	int size_x=img->width;
	for (int y = 0; y < img->height; y++)
	{
		for (int x = 0; x < img->width ; x++)
		{
			// pointer on the pixel of the images
			p = cvPtr2D (img, y, x, NULL);
			mat[x+y*size_x]=((double)*p-min_val)/max_val;
		}
	}
}

void convert_double_to_char(IplImage* img, double * mat,double max_val)
{
	int size_x=img->width;
	int size_y=img->height;
	uchar *p;
	double min_val_mat, max_val_mat;
	min_max_val(mat,min_val_mat,max_val_mat,size_x,size_y);
	for (int y = 0; y < img->height; y++)
	{
		for (int x = 0; x < img->width ; x++)
		{
			// pointer on the pixel of the images
			p = cvPtr2D (img, y, x, NULL);
			*p=((mat[x+y*size_x]-min_val_mat)/((max_val_mat-min_val_mat)))*max_val;
		}
	}
}

void smooth_img(double * img, double sigma, int size_x, int size_y)
{
	double s = 1.0/5.0;
	double t = sigma*sigma/2;
	int n_max = ceil(t/s);
	double * smooth= new double[size_y*size_x];
	for (int n=1 ; n < n_max ; n++){
		for (int y=1 ; y < size_y-1 ; y++){
			for (int x =1 ; x < size_x-1 ; x++){
				smooth[x+y*size_x]=img[x+y*size_x]+s*(img[x+y*size_x+1]+img[x+y*size_x-1]+img[x+y*size_x+size_x]+img[x+y*size_x-size_x]-4*img[x+y*size_x]);
			}
		}
		// borders
		for (int x =1 ; x < size_x-1 ; x++){
				smooth[x]=img[x]+s*(img[x+1]+img[x-1]+img[x+size_x]-3*img[x]);
				smooth[x+(size_y-1)*size_x]=img[x+(size_y-1)*size_x]+s*(img[x+(size_y-1)*size_x+1]+img[x+(size_y-1)*size_x-1]+img[x+(size_y-1)*size_x-size_x]-3*img[x+(size_y-1)*size_x]);
		}
		for (int y =1 ; y < size_y-1 ; y++){
				smooth[y*size_x]=img[y*size_x]+s*(img[y*size_x+1]+img[y*size_x+size_x]+img[y*size_x-size_x]-3*img[y*size_x]);
				smooth[y*size_x+(size_x-1)]=img[y*size_x+(size_x-1)]+s*(img[y*size_x+(size_x-1)-1]+img[(size_x-1)+y*size_x+size_x]+img[(size_x-1)+y*size_x-size_x]-3*img[(size_x-1)+y*size_x]);
		}
		//angles
		smooth[0]=img[0]+s*(img[1]+img[size_x]-2*img[0]);
		smooth[size_x-1]=img[size_x-1]+s*(img[size_x-2]+img[2*size_x-1]-2*img[size_x-1]);
		smooth[(size_y-1)*size_x]=img[(size_y-1)*size_x]+s*(img[(size_y-1)*size_x+1]+img[(size_y-2)*size_x]-2*img[(size_y-1)*size_x]);
		smooth[size_x*size_y-1]=img[size_x*size_y-1]+s*(img[size_x*size_y-2]+img[size_x*size_y-1-size_x]-2*img[size_x*size_y-1]);

		// copy
		for (int y=0 ; y < size_y ; y++){
			for (int x =0 ; x < size_x ; x++){
				img[x+y*size_x]=smooth[x+y*size_x];
			}
		}
	}
	// Normalization
	if (t!=0){
		for (int y=0 ; y < size_y ; y++){
			for (int x =0 ; x < size_x ; x++){
				img[x+y*size_x]=2*sqrt(M_PI*t)*img[x+y*size_x];
			}
		}
	}

	delete [] smooth;

}
