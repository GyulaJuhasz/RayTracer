#include <iostream>

#include "RealRoots.h"

using raytracer::RealRoots;

RealRoots::RealRoots() {
	numOfRoots = 0;
}

void RealRoots::AddRoot(double x0) {
	if (numOfRoots < MAX_POLINOM_DEGREE) {
		realRoots[numOfRoots++] = x0;
	}
}

int RealRoots::GetNumberOfRoots() { return numOfRoots; }

double RealRoots::GetRoot(int index) {
	if (index < 0 || index > numOfRoots) {
		throw "Wrong INDEX!";
	}

	return realRoots[index];
}

double RealRoots::GetMinPositiveRootIfAny() {
	if (numOfRoots == 0) return -1.0;

	bool found = false;
	double minPos = 99999999.0;
	double actualRoot;

	for (int i = 0; i < numOfRoots; i++) {
		actualRoot = realRoots[i];
		if (actualRoot > 0.0) {
			found = true;
		}

		if (actualRoot < minPos) {
			minPos = actualRoot;
		}
	}

	return (found ? minPos : -1.0);

}

RealRoots& RealRoots::operator=(const RealRoots& r) {
	if (this != &r) {
		numOfRoots = r.numOfRoots;
		for (int i = 0; i <= numOfRoots; i++) {
			realRoots[i] = r.realRoots[i];
		}
	}

	return *this;
}

std::ostream& raytracer::operator<<(std::ostream& os, const RealRoots& r) {
	if (r.numOfRoots == 0) {
		os << "NO REAL ROOTS!" << std::endl;
	}
	else {
		for (int i = 0; i < r.numOfRoots; i++) {
			os << "x" << (i + 1) << " = " << (r.realRoots[i]) << std::endl;
		}
	}

	return os;
}