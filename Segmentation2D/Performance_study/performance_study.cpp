/*
 * performance_study.cpp
 *
 *  Created on: 28 juin 2016
 *      Author: bastien.kovac
 */

#include "performance_study.h"

std::string get_line(std::string file_name, int nb_line) {

	std::string line;

	std::ifstream file(file_name.c_str());
	for (int i = 0 ; i < nb_line ; ++i) {
		getline(file, line);
	}
	getline(file, line);
	file.close();
	return line;

}

void edit_line(std::string file_name, std::string new_value, int nb_line) {

	std::vector<std::string> lines;
	std::string input;
	std::ifstream ifile(file_name.c_str());
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
	std::ofstream ofile(file_name.c_str());
	for (int i = 0 ; i < lines.size() ; i++) {
		ofile << lines[i] << std::endl;
	}

	ofile.close();

}

std::vector<double> test_line(double start, double end, double step,
		int argc, char ** argv, int nb_iterations, int line_to_modify) {

	int inc = start;
	int launch, stop;
	std::vector<double> results;

	if (inc >= 0) {

		do {
			edit_line("../Parameters.txt", std::to_string(inc), line_to_modify);
			launch = clock();
			main_logic(argc, argv, nb_iterations);
			stop = clock();
			results.push_back((double) (stop - launch) / double(CLOCKS_PER_SEC));
			std::cout << "\t" << inc << "/" << end << " done" << std::endl;
			inc += step;
		} while (!(inc > end));

	}

	return results;

}

void performance_test(int argc, char ** argv) {

	int starting_nb_ellipses, ending_nb_ellipses, step_nb_ellipses;

	int nb_iterations_start, nb_iterations_end, nb_iterations_step, nb_iterations;

	// we reset total_fit content
	std::ofstream total_fit_file_reset;
	total_fit_file_reset.open("../Total_Fit.txt");
	total_fit_file_reset.close();

	std::string st;
	std::ifstream ifile("../StudyParameters.txt", std::ios::in); // To open the file
	if (ifile) {
		ifile >> st;	ifile >> starting_nb_ellipses;
		ifile >> st;	ifile >> ending_nb_ellipses;
		ifile >> st;	ifile >> step_nb_ellipses;
		ifile >> st;	ifile >> nb_iterations_start;
		ifile >> st;	ifile >> nb_iterations_end;
		ifile >> st;	ifile >> nb_iterations_step;
	} else {
		std::cerr << "Error to open the file " << "StudyParameters.txt" << std::endl;
	}

	for (int nb_iterations = nb_iterations_start;
			nb_iterations <= nb_iterations_end; nb_iterations +=
					nb_iterations_step) {

		// Copy original file
		std::ifstream src("../Parameters_Model.txt", std::ios::binary);
		std::ofstream dst("../Parameters.txt", std::ios::binary);

		dst << src.rdbuf();

		src.close();
		dst.close();

		std::cout << "For " << nb_iterations << " iterations..." << std::endl;

		// For gradient
		edit_line("../Parameters.txt", "2", 12);
		// For monothread
		edit_line("../Parameters.txt", "1", 24);

		std::cout << "Testing number of ellipses for monothread ..." << std::endl;

		std::vector<double> grad_mono = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations,
				2);

		std::cout << "Done" << std::endl;

		// For multithread
		edit_line("../Parameters.txt", "0", 24);

		std::cout << "Testing number of ellipses for multithread ..." << std::endl;

		std::vector<double> grad_multi = test_line(starting_nb_ellipses,
				ending_nb_ellipses, step_nb_ellipses, argc, argv, nb_iterations,
				2);

		std::cout << "Done" << std::endl;

		std::vector<int> nb_ell;

		int inc = starting_nb_ellipses;

		do {
			nb_ell.push_back(inc);
			inc += step_nb_ellipses;
		} while (!(inc > ending_nb_ellipses));

		std::vector<std::string> lines;
		std::string line;

		std::string total_fit_file = "Total_Fit";

		std::vector<int> nb_it;
		std::vector<double> final_energy;

		for (int i = 1000; i <= 10000; i += 1000) {
			nb_it.push_back(i);
		}

		for (int i = 0; i < nb_ell.size(); i++) {
			line = (std::to_string(nb_ell[i]) + ";" + std::to_string(grad_mono[i]) + ";"
					+ get_line("../Total_Fit.txt", i) + ";"
					+ std::to_string(grad_multi[i]) + ";"
					+ get_line("../Total_Fit.txt", i + nb_ell.size()));
			lines.push_back(line);
		}

		//Write csv file
		std::ofstream csv_file;
		csv_file.open("../output.csv");

		csv_file << "Results for " << nb_iterations << " iterations" << std::endl;
		csv_file << std::endl;
		csv_file << ";gradient_mono;;gradient_multi" << std::endl;
		csv_file << std::endl;
		csv_file << "nb_ellipse;Seconds;Final energy;Seconds;Final energy"
				<< std::endl;
		for (int i = 0; i < lines.size(); i++) {
			csv_file << lines[i] << std::endl;
		}

		csv_file.close();

	}

}
