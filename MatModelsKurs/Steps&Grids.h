#ifndef STEPS_AND_GRIDS_H 
#define STEPS_AND_GRIDS_H 

#include <vector>

std::vector<double> getMainSteps(double lowerBound, double upperBound, int partsAmount);

std::vector<double> getMainGrid(double lowerBound, std::vector<double> mainSteps);

std::vector<double> getAuxiliarySteps(double lowerBound, std::vector<double> mainSteps);

std::vector<double> getAuxiliaryGrid(double lowerBound, std::vector<double> mainGrid);

#endif