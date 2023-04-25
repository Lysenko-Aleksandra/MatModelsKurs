#include "Steps&Grids.h"

std::vector<double> getMainSteps(double lowerBound, double upperBound, int partsAmount) {
	double h_equal = (upperBound - lowerBound) / partsAmount;//шаг основной сетки, сетка равномерная
	std::vector<double> mainSteps(partsAmount + 1, h_equal);///массив с шагами основной сетки

	return mainSteps;
}


std::vector<double> getMainGrid(double lowerBound, std::vector<double> mainSteps) {
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
			auxiliarySteps[i] = lowerBound + mainSteps[i + 1] / 2;
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

std::vector<double> getAuxiliaryGrid(double lowerBound, std::vector<double> mainGrid) {
	double gridSize = mainGrid.size() + 1;
	std::vector<double>auxGrid(gridSize);
	auxGrid[0] = lowerBound;
	auxGrid[1] = auxGrid[0] + (mainGrid[1] - mainGrid[0]) / 2;///
	for (int i = 2; i < gridSize-1; i++) {
		auxGrid[i] = auxGrid[i - 1] + (mainGrid[i] - mainGrid[i - 1]);
	}
	auxGrid[mainGrid.size()] = auxGrid[mainGrid.size() - 1] + (mainGrid[mainGrid.size() - 1] - mainGrid[mainGrid.size() - 2]) / 2;
	return auxGrid;
}