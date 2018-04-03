#include <iostream>

#include "MathUtils.h"
#include "Polinom.h"

using raytracer::Polinom;

Polinom::Polinom(int newDegree) {
	if (newDegree < 0 || newDegree > MAX_POLINOM_DEGREE) {
		throw "Invalid degree!";
	}

	degree = newDegree;
}

Polinom::Polinom(const Polinom& p) {
	degree = p.degree;
	for (int i = 0; i <= degree; i++) {
		coefficients[i] = p.coefficients[i];
	}
}

int Polinom::GetDegree() { return degree; }

double Polinom::GetCoefficient(int whichDegree) const {
	if (whichDegree < 0 || whichDegree > degree) {
		throw "Invalid degree!";
	}

	return coefficients[whichDegree];
}

void Polinom::SetCoeff(int whichDegree, double coeff) {
	if (whichDegree < 0 || whichDegree > degree) {
		throw "Invalid degree!";
	}

	coefficients[whichDegree] = coeff;
}

Polinom& Polinom::SetCoefficient(int whichDegree, double coeff) {
	if (whichDegree < 0 || whichDegree > degree) {
		throw "Invalid degree!";
	}

	coefficients[whichDegree] = coeff;
	return *this;
}

Polinom Polinom::divideByRoot(double x0) {
	if (degree == 0) {
		throw "Cannot divide 0 degree!";
	}

	Polinom result(degree - 1);

	for (int i = degree - 1; i >= 0; i--) {
		double ai = 0.0;

		for (int j = 0; j <= degree - 1 - i; j++) {
			double power;
			power = 1.0;

			for (int k = 0; k < (degree - 1 - i - j); k++) power *= x0;

			ai += GetCoefficient(degree - j) * power;
		}

		result.SetCoefficient(i, ai);
	}

	return result;
}

Polinom Polinom::GetDerived() {
	if (degree == 0) {
		Polinom result(0);
		result.SetCoefficient(0, 0.0);
		return result;
	}

	Polinom result(degree - 1);

	for (int i = degree - 1; i >= 0; i--) {
		result.SetCoefficient(i, GetCoefficient(i + 1) * (i + 1));
	}

	return result;
}

Polinom& Polinom::operator=(const Polinom& p) {
	if (this != &p) {
		degree = p.degree;
		for (int i = 0; i <= degree; i++) {
			coefficients[i] = p.coefficients[i];
		}
	}

	return *this;
}

double Polinom::GetValueAt(double x0) {
	double result = 0.0;
	double power = 1.0;

	for (int i = 0; i <= degree; i++) {
		if (i != 0) power *= x0;
		result += coefficients[i] * power;
	}

	return result;
}

raytracer::RealRoots Polinom::getRootsWithNewton(double guess) {
	RealRoots result;

	double prev = guess, fPrev, fDerivPrev;
	double actual;
	int goodTries = 0;
	int badTries = 0;

	double previousDistance = 9999999.0;
	double distance;

	Polinom basePoli(*this);
	Polinom derivedPoli = GetDerived();

	int maxRoots = basePoli.GetDegree();

	do {
		fPrev = basePoli.GetValueAt(prev);
		fDerivPrev = derivedPoli.GetValueAt(prev);

		if (raytracer::math::IsZero(fDerivPrev)) {
			return result;
		}

		actual = prev - (fPrev / fDerivPrev);

		distance = fabs(actual - prev);

		if (distance < numeric_epsilon) {
			goodTries++;
		}
		else if (distance > previousDistance) {
			badTries++;
		}

		previousDistance = distance;
		prev = actual;

		if (badTries > max_wrong_tries) {
			return result;
		}

		if (goodTries > min_good_tries) {
			result.AddRoot(actual);
			basePoli = basePoli.divideByRoot(actual);
			derivedPoli = basePoli.GetDerived();
			goodTries = 0;
			badTries = 0;
			prev = guess;
			previousDistance = 9999999.0;
		}

	} while (result.GetNumberOfRoots() < maxRoots);

	return result;

}

raytracer::RealRoots Polinom::getRootsWithBisection(Polinom polinom, double initA, double initB, raytracer::RealRoots *foundRoots) {
	RealRoots result;

	if (foundRoots != NULL) {
		result = *foundRoots;
	}

	int numOfSteps = 1000;

	double A = initA;
	double B = initB;
	double C;
	double step = (B - A) / numOfSteps;

	Polinom basePoli = polinom;

	int sigA = raytracer::math::Sig(basePoli.GetValueAt(A));
	int sigB = raytracer::math::Sig(basePoli.GetValueAt(B));
	int sigC;

	if (sigA == 0) {
		result.AddRoot(A);
		result = getRootsWithBisection(basePoli.divideByRoot(A), initA, initB, &result);
		return result;
	}

	if (sigB == 0) {
		result.AddRoot(B);
		result = getRootsWithBisection(basePoli.divideByRoot(B), initA, initB, &result);
		return result;
	}

	while (sigA == sigB && A < B) {
		A += step;
		sigA = raytracer::math::Sig(basePoli.GetValueAt(A));
	}

	if (!(A < B)) return result;

	// There is at least one root

	double root;
	do {
		if (sigA == 0) {
			root = A;
			break;
		}

		if (sigB == 0) {
			root = B;
			break;
		}

		C = 0.5 * (A + B);
		sigC = raytracer::math::Sig(basePoli.GetValueAt(C));

		if (sigC == 0) {
			root = C;
			break;
		}

		if (abs(A - B) < 0.00001) {
			root = C;
			break;
		}

		if (sigC == sigA) {
			A = C;
			sigA = sigC;
		}
		else if (sigC == sigB) {
			B = C;
			sigB = sigC;
		}

	} while (true);

	result.AddRoot(root);
	result = getRootsWithBisection(basePoli.divideByRoot(root), initA, initB, &result);
	return result;
}

const double Polinom::numeric_epsilon = 0.000001;
const int Polinom::max_wrong_tries = 10;
const int Polinom::min_good_tries = 3;

std::ostream& raytracer::operator<<(std::ostream& os, const Polinom& p) {
	for (int i = p.degree; i >= 0; i--) {
		double next;
		os << (i == p.degree ? p.GetCoefficient(i) : fabs(p.GetCoefficient(i)));
		if (i > 0) {
			next = p.GetCoefficient(i - 1);
			os << " * x";
			if (i > 1) {
				os << "^" << i;
			}
			os << (next > 0.0 ? " + " : " - ");
		}
	}
	return os;
}
