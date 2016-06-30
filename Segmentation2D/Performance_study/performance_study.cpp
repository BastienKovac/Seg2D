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

string get_line(string file_name, int nb_line) {

	string line;

	ifstream file(file_name.c_str());
	for (int i = 0 ; i < nb_line ; ++i) {
		getline(file, line);
	}
	getline(file, line);
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

	int inc = start;
	int launch, stop;
	vector<double> results;

	if (inc >= 0) {

		do {
			edit_line("../Parameters.txt", to_string(inc), line_to_modify);
			launch = clock();
			main_logic(argc, argv, nb_iterations);
			stop = clock();
			results.push_back((double) (stop - launch) / double(CLOCKS_PER_SEC));
			cout << "\t" << inc << "/" << end << " done" << endl;
			inc += step;
		} while (!(inc > end));

	}

	return results;

}

void performance_test(int argc, char ** argv) {

	int starting_nb_ellipses, ending_nb_ellipses, step_nb_ellipses;

	int nb_iterations_start, nb_iterations_end, nb_iterations_step, nb_iterations;

	// we reset total_fit content
	ofstream total_fit_file_reset;
	total_fit_file_reset.open("../Total_Fit.txt");
	total_fit_file_reset.close();

	string st;
	ifstream ifile("../StudyParameters.txt", ios::in); // To open the file
	if (ifile) {
		ifile >> st;	ifile >> starting_nb_ellipses;
		ifile >> st;	ifile >> ending_nb_ellipses;
		ifile >> st;	ifile >> step_nb_ellipses;
		ifile >> st;	ifile >> nb_iterations_start;
		ifile >> st;	ifile >> nb_iterations_end;
		ifile >> st;	ifile >> nb_iterations_step;
	} else {
		cerr << "Error to open the file " << "StudyParameters.txt" << endl;
	}

	for (int nb_iterations = nb_iterations_start;
			nb_iterations <= nb_iterations_end; nb_iterations +=
					nb_iterations_step) {

		// Copy original file
		ifstream src("../Parameters_Model.txt", ios::binary);
		ofstream dst("../Parameters.txt", ios::binary);

		dst << src.rdbuf();

		src.close();
		dst.close();

		cout << "For " << nb_iterations << " iterations..." << endl;

		// For gradient
		edit_line("../Parameters.txt", "2", 12);
		// For monothread
		edit_line("../Parameters.txt", "1", 24);

		cout << "Testing number of ellipses for monothread ..." << endl;

		vector<double> grad_mono = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations,
				2);

		cout << "Done" << endl;

		// For multithread
		edit_line("../Parameters.txt", "0", 24);

		cout << "Testing number of ellipses for multithread ..." << endl;

		vector<double> grad_multi = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations,
				2);

		cout << "Done" << endl;

		vector<int> nb_ell;

		int inc = starting_nb_ellipses;

		do {
			nb_ell.push_back(inc);
			inc += step_nb_ellipses;
		} while (!(inc > ending_nb_ellipses));

		vector<string> lines;
		string line;

		string total_fit_file = "Total_Fit";

		vector<int> nb_it;
		vector<double> final_energy;

		for (int i = 1000; i <= 10000; i += 1000) {
			nb_it.push_back(i);
		}

		for (int i = 0; i < nb_ell.size(); i++) {
			line = (to_string(nb_ell[i]) + ";" + to_string(grad_mono[i]) + ";"
					+ get_line("../Total_Fit.txt", i) + ";"
					+ to_string(grad_multi[i]) + ";"
					+ get_line("../Total_Fit.txt", i + nb_ell.size()));
			lines.push_back(line);
		}

		//Write csv file
		ofstream csv_file;
		csv_file.open("../output.csv");

		csv_file << "Results for " << nb_iterations << " iterations" << endl;
		csv_file << endl;
		csv_file << ";gradient_mono;;gradient_multi" << endl;
		csv_file << endl;
		csv_file << "nb_ellipse;Seconds;Final energy;Seconds;Final energy"
				<< endl;
		for (int i = 0; i < lines.size(); i++) {
			csv_file << lines[i] << endl;
		}

		csv_file.close();

	}

}
