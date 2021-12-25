#ifndef DOUBLE_COMPARATOR_H
#define DOUBLE_COMPARATOR_H

#include <cmath>

namespace raytracer {

	namespace math {

		const double epsilon = 0.000001;

		inline bool Equals(double a, double b) {
			return fabs(a - b) < epsilon;
		}

		inline bool NotEquals(double a, double b) {
			return !Equals(a, b);
		}

		inline bool IsZero(double a) {
			return Equals(a, 0.0);
		}

		inline int Sig(double a) {
			if (a > epsilon / 10) {
				return 1;
			}
			else if (a < -epsilon / 10) {
				return -1;
			}
			else {
				return 0;
			}
		}
	}

}

#endif // !DOUBLE_COMPARATOR_H