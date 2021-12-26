#ifndef CONE_H
#define CONE_H

#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Cone : public raytracer::Object {
		raytracer::Vector baseCenter;
		double baseRadius;
		double height;
		bool isSolid;

		// Calculated fields
		raytracer::Vector vertex;
		double baseRadiusPerHeight;
		double sinAlpha;
		double cosAlpha;
		double oneMinusCosAlpha;

		void CalculateValues();
		raytracer::Vector GetSurfaceNormal(raytracer::Vector surfacePoint);
	public:
		Cone();

		Cone& BaseCenter(raytracer::Vector newBaseCenter);

		Cone& BaseRadius(double newBaseRadius);

		Cone& Height(double newBaseHeight);

		Cone& IsSolid(bool newSolid);

	protected:
		Intersection SpecificIntersect(raytracer::Ray ray);
	};
};

#endif // !CONE_H