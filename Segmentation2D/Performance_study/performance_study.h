/*
 * performance_study.h
 *
 *  Created on: 28 juin 2016
 *      Author: bastien.kovac
 */

#include <string>
#undef max
#include <vector>

#ifndef PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_
#define PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_

// Loops the main algorithm according to the StudyParameters file and print the results in an output csv file
void performance_test(int argc, char ** argv);

void edit_line(std::string file_name, std::string new_value, int nb_line);

std::string to_string(double nb);

std::vector<double> test_line(double start, double end, double step,
		int argc, char ** argv, int nb_iterations, int line_to_modify);

#endif /* PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_ */
