/*
 * performance_study.cpp
 *
 *  Created on: 28 juin 2016
 *      Author: bastien.kovac
 */

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#undef max
#include <vector>
#include "../Main_prog/Segmentation_prog.h"
#include "performance_study.h"

using namespace std;

string to_string(double nb) {
	return static_cast<ostringstream*>( &(ostringstream() << nb) )->str();
}

string getline(string file_name, int nb_line) {

	string line = "";
	ifstream file(file_name.c_str());
	int nbline = 0;

	while(!file.eof()) {
		getline(file, line);
		if (nbline == nb_line - 1) {
			break;
		}
		nbline++;
	}

	file.close();

	return line;

}

void edit_line(std::string file_name, std::string new_value, int nb_line) {

	vector<string> lines;
	string input;
	ifstream ifile(file_name.c_str());
	while (getline(ifile, input)) {
		lines.push_back(input);
	}

	for (int i = 0 ; i < lines.size() ; i++) {
		if (i == nb_line - 1) {
			lines[i] = new_value;
			break;
		}
	}

	ifile.close();
	ofstream ofile(file_name.c_str());
	for (int i = 0 ; i < lines.size() ; i++) {
		ofile << lines[i] << endl;
	}

	ofile.close();

}

std::vector<double> test_line(double start, double end, double step,
		int argc, char ** argv, int nb_iterations, int line_to_modify) {

	// Copy original file
	ifstream src("Parameters_Model.txt", ios::binary);
	ofstream dst("Parameters.txt", ios::binary);

	dst << src.rdbuf();

	src.close();
	dst.close();

	int inc = start;

	int launch, stop;

	vector<double> results;

	if (inc >= 0) {

		do {
			edit_line("Parameters.txt", to_string(inc), line_to_modify);
			launch = clock();
			main_logic(argc, argv, nb_iterations);
			stop = clock();
			inc += step;
			results.push_back((double) (stop - launch) / double(CLOCKS_PER_SEC));
		} while (!(inc > end));

	}

	return results;

}

void performance_test(int argc, char ** argv) {

	int starting_nb_ellipses, ending_nb_ellipses, step_nb_ellipses;

	int nb_iterations;

	// we reset total_fit content
	ofstream total_fit_file_reset;
	total_fit_file_reset.open("Total_Fit.txt");
	total_fit_file_reset.close();

	string st;
	ifstream ifile("StudyParameters.txt", ios::in); // To open the file
	if (ifile) {
		ifile >> st;	ifile >> starting_nb_ellipses;
		ifile >> st;	ifile >> ending_nb_ellipses;
		ifile >> st;	ifile >> step_nb_ellipses;
		ifile >> st;	ifile >> nb_iterations;
	} else {
		cerr << "Error to open the file " << "StudyParameters.txt" << endl;
	}

	// For first method
	edit_line("Parameters.txt", "1", 12);
	// For monothread
	edit_line("Parameters.txt", "1", 24);

	vector<double> contrast_mono = test_line(starting_nb_ellipses,
			ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations, 2);

	// For multithread
	edit_line("Parameters.txt", "0", 24);

	vector<double> contrast_multi = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations, 2);

	// For second method
	edit_line("Parameters.txt", "2", 12);
	// For monothread
	edit_line("Parameters.txt", "1", 24);

	vector<double> grad_mono = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations, 2);

	// For multithread
	edit_line("Parameters.txt", "0", 24);

	vector<double> grad_multi = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations, 2);

	// Copy original file
	ifstream src("Parameters_Model.txt", ios::binary);
	ofstream dst("Parameters.txt", ios::binary);

	dst << src.rdbuf();

	vector<int> nb_ell;

	int inc = starting_nb_ellipses;

	do {
		nb_ell.push_back(inc);
		inc += step_nb_ellipses;
	} while (!(inc > ending_nb_ellipses));

	vector<string> lines;
	string line;

	string total_fit_file = "Total_Fit";

	for (int i = 0 ; i < nb_ell.size() ; i++) {
		line = (to_string(nb_ell[i]) + ";"
				+ to_string(contrast_mono[i]) + ";"
				+ getline(total_fit_file, i) + ";"
				+ to_string(contrast_multi[i]) + ";"
				+ getline(total_fit_file, i + 0 * nb_ell.size()) + ";"
				+ to_string(grad_mono[i]) + ";"
				+ getline(total_fit_file, i + 2 * nb_ell.size()) + ";"
				+ to_string(grad_multi[i])
				+ getline(total_fit_file, i + 3 * nb_ell.size()) + ";");
		lines.push_back(line);
	}

	nb_ell.clear();
	vector<int>().swap(nb_ell);

	contrast_mono.clear();
	vector<double>().swap(contrast_mono);

	contrast_multi.clear();
	vector<double>().swap(contrast_multi);

	grad_mono.clear();
	vector<double>().swap(grad_mono);

	grad_multi.clear();
	vector<double>().swap(grad_multi);

	//Write csv file
	ofstream csv_file;
	csv_file.open("output.csv");

	csv_file << "Results" << endl;
	csv_file << endl;
	csv_file << ";contrast_mono;;contrast_multi;;gradient_mono;;gradient_multi" << endl;
	csv_file << endl;
	csv_file << "nb_ellipse; seconds ; energy ; seconds ; energy ; seconds ; energy ; seconds ; energy" << endl;
	for (int i = 0 ; i < lines.size() ; i++) {
		csv_file << lines[i] << endl;
	}

	csv_file.close();

}
