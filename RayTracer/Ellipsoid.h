#ifndef ELLIPSIOD_H
#define ELLIPSIOD_H

#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Ellipsoid : public raytracer::Object {
		raytracer::Vector center;
		double _A, _B, _C;
	public:
		Ellipsoid();

		Ellipsoid& Center(raytracer::Vector newCenter);

		raytracer::Vector GetCenter();

		Ellipsoid& A(double newA);

		double GetA();

		Ellipsoid& B(double newB);

		double GetB();

		Ellipsoid& C(double newC);

		double GetC();

		Intersection Intersect(raytracer::Ray ray);
	};

}

#endif // !ELLIPSIOD_H