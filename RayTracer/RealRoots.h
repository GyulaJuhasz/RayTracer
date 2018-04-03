#ifndef REAL_ROOTS_H
#define REAL_ROOTS_H

#include <iosfwd>

#include "Config.h"

namespace raytracer {

	class RealRoots {
		double realRoots[MAX_POLINOM_DEGREE];
		int numOfRoots;
	public:
		RealRoots();

		void AddRoot(double x0);

		int GetNumberOfRoots();

		double GetRoot(int index);

		double GetMinPositiveRootIfAny();

		RealRoots& operator=(const RealRoots& r);

		friend std::ostream& operator<<(std::ostream& os, const RealRoots& r);

	};

	std::ostream& operator<<(std::ostream& os, const RealRoots& r);

}

#endif // !REAL_ROOTS_H
