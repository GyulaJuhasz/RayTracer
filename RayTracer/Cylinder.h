#ifndef CYLINDER_H
#define CYLINDER_H

#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Cylinder : public raytracer::Object {
		raytracer::Vector basePoint;
		double radius;
		double height;
		bool solid;
	public:
		Cylinder();

		Cylinder& BasePoint(raytracer::Vector newBasePoint);

		raytracer::Vector GetBasePoint();

		Cylinder& Radius(double newRadius);

		double GetRadius();

		Cylinder& Height(double newHeight);

		double GetHeight();

		Cylinder& Solid(bool newSolid = true);

		bool IsSolid();

	protected:
		Intersection SpecificIntersect(raytracer::Ray ray);
	};

}

#endif // !CYLINDER_H