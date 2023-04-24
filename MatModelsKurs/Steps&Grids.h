#ifndef STEPS_AND_GRIDS_H 
#define STEPS_AND_GRIDS_H 

#include <vector>

std::vector<double> getMainSteps(double lowerBound, double upperBound, int partsAmount) {
	double h_equal = (upperBound - lowerBound) / partsAmount;//шаг основной сетки, сетка равномерная
	std::vector<double> mainSteps(partsAmount + 1, h_equal);///массив с шагами основной сетки

	return mainSteps;
}


std::vector<double> getMainGrid(double lowerBound,std::vector<double> mainSteps) {
	std::vector<double> mainGrid(mainSteps.size());///по основной сетке(допускается случай с неоднородной сеткой)
	for (int i = 0; i < mainSteps.size(); i++) {
		if (i == 0) {
			mainGrid[i] = lowerBound;
		}
		else {
			mainGrid[i] = mainGrid[i - 1] + mainSteps[i];
		}
	}

	return mainGrid;
}

std::vector<double> getAuxiliarySteps(double lowerBound, std::vector<double> mainSteps) {
	int partsAmount = mainSteps.size();
	std::vector<double> auxiliarySteps(partsAmount);//массив с шагами дополнительной сетки
	for (int i = 0; i < partsAmount; i++) {
		if (i == 0) {
			auxiliarySteps[i] = lowerBound+mainSteps[i + 1] / 2;
		}
		else if (i == partsAmount) {
			auxiliarySteps[i] = mainSteps[i] / 2;
		}
		else {
			auxiliarySteps[i] = (mainSteps[i] + mainSteps[i + 1]) / 2;
		}
	}
	return auxiliarySteps;
}

std::vector<double> getAuxiliaryMinusHalfGrid(double lowerBound, std::vector<double> mainGrid) {
	int partsAmount = mainGrid.size();
	std::vector<double> auxiliaryMinusHalfGrid(partsAmount);///r -1/2
	for (int i = 1; i <partsAmount; i++) {
		auxiliaryMinusHalfGrid[i] = (mainGrid[i] + mainGrid[i - 1]) / 2;
	}
	return auxiliaryMinusHalfGrid;
}

std::vector<double> getAuxiliaryPlusHalfGrid(double lowerBound, std::vector<double> mainGrid) {
	int partsAmount = mainGrid.size();
	std::vector<double> auxiliaryPlusHalfGrid(partsAmount);///r +1/2
	for (int i = 0; i < partsAmount-1; i++) {
		auxiliaryPlusHalfGrid[i] = (mainGrid[i] + mainGrid[i + 1]) / 2;
	}
	return auxiliaryPlusHalfGrid;
}

#endif