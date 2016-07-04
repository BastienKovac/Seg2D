/*
 * performance_study.h
 *
 *  Created on: 28 juin 2016
 *      Author: bastien.kovac
 */


#ifndef PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_
#define PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>
#include <vector>
#include "../Main_prog/Segmentation_prog.h"

// Loops the main algorithm according to the StudyParameters file and print the results in an output csv file
void performance_test(int argc, char ** argv);

void edit_line(std::string file_name, std::string new_value, int nb_line);

std::string get_line(std::string file_name, int nb_line);

std::vector<double> test_line(double start, double end, double step,
		int argc, char ** argv, int nb_iterations, int line_to_modify);

#endif /* PERFORMANCE_STUDY_PERFORMANCE_STUDY_H_ */
