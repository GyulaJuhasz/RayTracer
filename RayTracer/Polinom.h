#ifndef POLINOM_H
#define POLINOM_H

#include <iosfwd>

#include "Config.h"
#include "RealRoots.h"

namespace raytracer {

	class Polinom {
		double coefficients[MAX_POLINOM_DEGREE];
		int degree;

		static const double numeric_epsilon;
		static const int max_wrong_tries;
		static const int min_good_tries;

	public:
		Polinom(int newDegree);

		Polinom(const Polinom& p);

		int GetDegree();

		double GetCoefficient(int whichDegree) const;

		void SetCoeff(int whichDegree, double coeff);

		Polinom& SetCoefficient(int whichDegree, double coeff);

		Polinom divideByRoot(double x0);

		Polinom GetDerived();

		Polinom& operator=(const Polinom& p);

		double GetValueAt(double x0);

		raytracer::RealRoots getRootsWithNewton(double guess);

		raytracer::RealRoots getRootsWithBisection(Polinom polinom, double initA, double initB, raytracer::RealRoots *foundRoots);

		friend std::ostream& operator<<(std::ostream& os, const Polinom& p);
	};

}

#endif // !POLINOM_H