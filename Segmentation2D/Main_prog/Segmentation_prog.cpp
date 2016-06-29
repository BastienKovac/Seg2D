/* ---------- Principal program ------------
 file Segmentation_prog.cpp  Modified on 09/07/1012
 source file

 Description:
 This program compute the segmentation of an image
 using the graphs cuts
 PARAMETERS
 - Name of the image
 RETURN
 - Segmentation_name.txt -> file which contains the parameters of the Ellips configuration
 - Segmentation_name.png -> image of the segmentation
 ---------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include "../Class_Ellips/Ellips.h"
#include "../Graph_Cut/graph.h"
#include "../Class_Configuration/Configuration.h"
#include "../other_functions/other_functions.h"
#include "../Performance_study/performance_study.h"
#include "Segmentation_prog.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <highgui.h>
#include <cv.h>
#include <cstdlib>

#include <blas.h>
#include <lapacke.h>
#include <omp.h>

using namespace std;

int using_multithread = 1;

int main_logic(int argc, char** argv, int nb_iterations = 0) {

	clock_t start, end;
	double cpuTime;
	time_t t = time(NULL);

	start = clock();

	// test the good number of arguments
	if (argc != 2 && argc != 3) {
		cout << "Wrong number of arguments" << endl;
		cout << "Usage: " << argv[0] << "  Image_name" << endl;
		return 1;
	}

	//////////////////////////////////////////////////////////////
	///////                 DEFINITIONS                   ////////
	//////////////////////////////////////////////////////////////

	// Definition of an infinity constant which is equal to the maximal float representable
	float inf = FLT_MAX;

	// Definition of a Graph type
	typedef Graph<float, float, float> GraphType;

	srand(time(NULL)); //  /!\ Warning /!\ don't forget this line into the main program to initialize rand !!!

	int nb_ell, nb_cell; // number of Ellipse per configuration // estimation of the number of cells
	int nb_ell_dont_accepted, nb_ell_tot, nb_ell_config, nb_ell_new_config;
	double a_min, a_max, min_val, max_val;
	float flow;
	int choice, i, j, stop, nb_threads;
	Configuration config, new_config, temp;
	double d, sigma; //acceptance threshold
	double epsilon; // For gradient method

	// Definition of the image of type IplImage
	IplImage* img = NULL;
	IplImage* print = NULL;
	const char* window_title = "Segmentation display";
	const char* window_title_gradx = "Gradient x";
	const char* window_title_grady = "Gradient y";
	const char* window_title_smooth = "Smoothed image";

	// Loading the image
	img = cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	if (img == NULL) {
		cout << "couldn't open image file: " << argv[1] << endl;
		return 1;
	}

	string st, name;
	// Reading of the parameters in the file Parameters.txt
	ifstream ifile("Parameters.txt", ios::in); // To open the file
	if (ifile) {
		ifile >> st;	ifile >> nb_ell;
		ifile >> st;	ifile >> nb_ell_dont_accepted;
		ifile >> st;	ifile >> a_min;
		ifile >> st;	ifile >> a_max;
		ifile >> st;	ifile >> nb_cell;
		ifile >> st;	ifile >> choice;
		ifile >> st;	ifile >> name;
		ifile >> st;	ifile >> d;
		ifile >> st;	ifile >> stop;
		ifile >> st;	ifile >> sigma;
		ifile >> st;	ifile >> epsilon;
		ifile >> st;	ifile >> nb_threads;
		ifile.close(); // To close the file
	} else {
		cerr << "Error to open the file " << "Parameters.txt" << endl;
		return 0;
	} // if(file)

	using_multithread = 1;
	if (nb_threads == 1) {
		using_multithread = 0;
	} else {
		if (nb_threads != 0) {
			omp_set_num_threads(nb_threads);
		}
	}

	// Size of the image
	int size_x = cvGetSize(img).width;
	int size_y = cvGetSize(img).height;

	// Definition of matrix to convert char in double
	double * im = new double[size_y * size_x];
	double * gradx = new double[size_y * size_x];
	double * grady = new double[size_y * size_x];

	convert_char_to_double(img, im);
	min_max_val(img, min_val, max_val);

	if (choice == 2) {
		// smooth the image
		smooth_img(im, sigma, size_x, size_y);

		if (argc == 2) {
			print = cvCreateImage(cvGetSize(img), 8, 1);
			convert_double_to_char(print, im, max_val);
			cvNamedWindow(window_title_smooth, CV_WINDOW_AUTOSIZE);
			cvShowImage(window_title_smooth, print);
			cout << "Press enter to continue" << endl;
			cvWaitKey(0);
			cvDestroyAllWindows();
			cvReleaseImage(&print);
		}

		// Gradient of the image
		grad_X(im, gradx, size_x, size_y);
		grad_Y(im, grady, size_x, size_y);

		if (argc == 2) {
			print = cvCreateImage(cvGetSize(img), 8, 1);
			convert_double_to_char(print, gradx, max_val);
			cvNamedWindow(window_title_gradx, CV_WINDOW_AUTOSIZE);
			cvShowImage(window_title_gradx, print);
			cout << "Press enter to continue" << endl;
			cvWaitKey(0);
			cvDestroyAllWindows();
			cvReleaseImage(&print);


			print = cvCreateImage(cvGetSize(img), 8, 1);
			convert_double_to_char(print, grady, max_val);
			cvNamedWindow(window_title_grady, CV_WINDOW_AUTOSIZE);
			cvShowImage(window_title_grady, print);
			cout << "Press enter to continue" << endl;
			cvWaitKey(0);
			cvDestroyAllWindows();
			cvReleaseImage(&print);
		}

	}

	// Definition of a pointer on a GraphType
	GraphType *g = new GraphType(/*estimated # of nodes*/4 * nb_ell, /*estimated # of edges*/
			2 * nb_ell * 4);

	// Generation of a first configuration
	if (choice == 1) {
		config = Configuration(a_min, a_max, size_x, size_y, nb_ell,
				nb_ell_dont_accepted, im, d); // Version 1
	}
	if (choice == 2) {
		config = Configuration(a_min, a_max, size_x, size_y, nb_ell,
				nb_ell_dont_accepted, gradx, grady, M_PI / 12, epsilon, d); // Version 2
	}

	try {
		// grid
		int nb_rows = ceil(size_y / a_max);
		int nb_col = ceil(size_x / a_max);

		// To stop the program
		bool accepted = false;
		int num_not_acc = 0;
		int k = 1;

<<<<<<< HEAD
		cvNamedWindow(window_title, CV_WINDOW_AUTOSIZE);
		// while(num_not_acc<stop)
		while (num_not_acc < stop) {
=======
		cvNamedWindow (window_title, CV_WINDOW_AUTOSIZE);

		while(num_not_acc<stop){
>>>>>>> branch 'master' of https://github.com/BastienKovac/Segmentation2D

			//---- Generation of a second configuration
			if (choice == 1) {
				new_config = Configuration(a_min, a_max, size_x, size_y, nb_ell,
						nb_ell_dont_accepted, im, d); // Version 1
			}
			if (choice == 2) {
				new_config = Configuration(a_min, a_max, size_x, size_y, nb_ell,
						nb_ell_dont_accepted, gradx, grady, M_PI / 12, epsilon,
						d); // Version 2
			}

			// number of Ellipses in config
			nb_ell_config = config.get_nb_Ellipses();
			// number of Ellipses in new_config
			nb_ell_new_config = new_config.get_nb_Ellipses();
			// total number of Ellipse (with the two configurations)
			nb_ell_tot = nb_ell_config + nb_ell_new_config;

			if (fmod(float(k), 1000) == 0) {
				cout << "Iteration : " << k << endl;
				cout << "Total execution time : " << difftime(time(NULL), t)
						<< "s" << endl;
				string file_name = "Segmentation_" + name + ".txt";

				config.save_config(file_name);
				print = cvCreateImage(cvGetSize(img), 8, 3);
				cvCvtColor(img, print, CV_GRAY2BGR);

				for (int z = 0; z < nb_ell_config; z++) {
					cvEllipse(print,
							cvPoint(config.get_Ellips(z).get_cx(),
									config.get_Ellips(z).get_cy()),
							cvSize(config.get_Ellips(z).get_a(),
									config.get_Ellips(z).get_b()),
							-config.get_Ellips(z).get_theta() * 360
									/ (2 * M_PI), 0, 360, CV_RGB(0, 0, 255), 1,
							8, 0);
				}
				cvShowImage(window_title, print);
				cvWaitKey(1);
				cvReleaseImage(&print);

			}

			//---- Add nodes to the graph
			g->add_node(nb_ell_tot);

			//---- Add the weights of the different edges
			for (i = 0; i < nb_ell_config; i++) {
				g->add_tweights(i, /* capacities */config.get_data_fit(i),
						1 - config.get_data_fit(i));
				for (j = 0; j < nb_ell_new_config; j++) {
					if (is_neighbor(config.get_position(i),
							new_config.get_position(j), nb_rows, nb_col,
							a_max)) {
						if (intersect(config.get_Ellips(i),
								new_config.get_Ellips(j))) {
							{
								g->add_edge(i, nb_ell_config + j, /* capacities */
										inf, 0);
							}
						}
					}
				}
			}

			#pragma omp parallel for if(using_multithread)
			for (i = 0; i < nb_ell_new_config; i++) {
				g->add_tweights(nb_ell_config + i, /* capacities */
						1 - new_config.get_data_fit(i),
						new_config.get_data_fit(i));
			}

			//---- Compute the max flow of the graph
			flow = g->maxflow();

			//---- keep the good Ellipses
			accepted = false;
			temp = config;
			i = 0;
			while ((g->what_segment(i) == GraphType::SINK) & (i < nb_ell_config)) {
				i++;
			}
			if (i < nb_ell_config) {
				config = Configuration(temp.get_Ellips(i), temp.get_position(i),
						temp.get_data_fit(i), nb_ell_tot);
				i++;
				for (; i < nb_ell_config; i++) {
					if (g->what_segment(i) == GraphType::SOURCE) {
						config.add_Ellips(temp.get_Ellips(i),
								temp.get_position(i), temp.get_data_fit(i));
					}
				}
				for (; i < nb_ell_tot; i++) {
					if (g->what_segment(i) == GraphType::SINK) {
						config.add_Ellips(
								new_config.get_Ellips(i - nb_ell_config),
								new_config.get_position(i - nb_ell_config),
								new_config.get_data_fit(i - nb_ell_config));
						accepted = true;
					}
				}
			} else {
				while ((i < nb_ell_tot)
						& (g->what_segment(i) == GraphType::SOURCE)) {
					i++;
					if (i == nb_ell_tot)
						break;
				}
				if (i < nb_ell_tot) {
					config = Configuration(
							new_config.get_Ellips(i - nb_ell_config),
							new_config.get_position(i - nb_ell_config),
							new_config.get_data_fit(i - nb_ell_config),
							nb_ell_tot);
					i++;
					accepted = true;
					for (; i < nb_ell_tot; i++) {
						if (g->what_segment(i) == GraphType::SINK) {
							config.add_Ellips(
									new_config.get_Ellips(i - nb_ell_config),
									new_config.get_position(i - nb_ell_config),
									new_config.get_data_fit(i - nb_ell_config));
						}
					}
				} else {
					// None of the node are kept so we generate a new_configuration
					if (choice == 1) {
						config = Configuration(a_min, a_max, size_x, size_y,
								nb_ell, nb_ell_dont_accepted, im, d); // Version 1
					}
					if (choice == 2) {
						config = Configuration(a_min, a_max, size_x, size_y,
								nb_ell, nb_ell_dont_accepted, gradx, grady,
								M_PI / 12, 0.001, d); // Version 2
					}

				} // if (i<nb_ell_tot)
			} // if (i<nb_ell_config)

			// reset the graph
			g->reset();

			if (accepted == true) {
				num_not_acc = 0;
			} else {
				num_not_acc++;
			}
			k++;

			if (nb_iterations != 0 && k == nb_iterations) {
				break;
			}

		} // while(num_not_acc<30000)

		if (argc == 3) {
			double total_fit = config.get_data_fit_total();
			ofstream total_fit_file;
			total_fit_file.open("Total_Fit.txt", ios_base::app);
			total_fit_file << total_fit << endl;
			total_fit_file.close();
		}

		// We save the segmentation
		string file_name = "Segmentation_" + name + ".txt";
		config.save_config(file_name);
		cout << "\nNumber of Ellipses found : " << config.get_nb_Ellipses()
				<< endl;

		// We save an image
		file_name = "Segmentation_" + name + ".png";
		print = cvCreateImage(cvGetSize(img), 8, 3);
		cvCvtColor(img, print, CV_GRAY2BGR);
		for (int z = 0; z < nb_ell_config; z++) {
			cvEllipse(print,
					cvPoint(config.get_Ellips(z).get_cx(),
							config.get_Ellips(z).get_cy()),
					cvSize(config.get_Ellips(z).get_a(),
							config.get_Ellips(z).get_b()),
					-config.get_Ellips(z).get_theta() * 360 / (2 * M_PI), 0,
					360, CV_RGB(0, 0, 255), 1, 8, 0);
		}
		cvSaveImage(file_name.c_str(), print);
		cvWaitKey(1);
		cvReleaseImage(&print);

		end = clock();
		// compute the executive time
		cpuTime = (end - start) / (CLOCKS_PER_SEC);
		double hours = (cpuTime / double(60)) / double(60);
		double minuts = (hours - floor(hours)) * double(60);
		double second = (minuts - floor(minuts)) * double(60);
		cout << "\nCPU Time : " << floor(hours) << " hours " << floor(minuts)
				<< " minuts " << second << " seconds " << endl;

		// Clean of the memory
		cvReleaseImage(&img);
		delete g;
	} catch (exception const& e) // We catch the exceptions.
	{
		cout << "ERREUR : " << e.what() << endl; // We plot the exceptions.
	}

	return 1;
}

int main(int argc, char ** argv) {

	if (argc == 2) {
		main_logic(argc, argv);
	} else {
		performance_test(argc, argv);
	}

	return 1;
}
