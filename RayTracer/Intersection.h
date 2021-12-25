#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "Vector.h"

namespace raytracer {

	class Intersection {
		raytracer::Vector intersectPoint;
		raytracer::Vector normal;
		double t;
		int index;

	public:
		Intersection();

		Intersection& Param(double param);

		double GetParam();

		Intersection& IntersectPoint(raytracer::Vector mp);

		raytracer::Vector GetIntersectPoint();

		Intersection& Normal(raytracer::Vector norm);

		raytracer::Vector GetNormal();

		Intersection& Index(int newIndex);

		int GetIndex();

		Intersection& operator=(const Intersection& other);

	};

}

#endif // !INTERSECTION_H