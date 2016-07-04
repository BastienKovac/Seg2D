/*
 * Segmentation_prog.h
 *
 *  Created on: 28 juin 2016
 *      Author: bastien.kovac
 */

#ifndef SEGMENTATION_PROG_H_
#define SEGMENTATION_PROG_H_

#include "../Class_Ellips/Ellips.h"
#include "../Graph_Cut/graph.h"
#include "../Class_Configuration/Configuration.h"
#include "../other_functions/other_functions.h"
#include "../Performance_study/performance_study.h"
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <float.h>
#include <cstdlib>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>

extern int using_multithread;

int main_logic(int argc, char** argv, int nb_iterations);


#endif /* SEGMENTATION_PROG_H_ */
