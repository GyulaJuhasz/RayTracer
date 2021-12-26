#ifndef HEART_SHAPE_H
#define HEART_SHAPE_H

#include "Cube.h"
#include "Intersection.h"
#include "Object.h"
#include "Ray.h"

namespace raytracer {

	class HeartShape : public raytracer::Object {

		raytracer::Cube cube;

		double minX;
		double maxX;
		double minY;
		double maxY;
		double minZ;
		double maxZ;

	public:
		HeartShape();

	protected:
		raytracer::Intersection SpecificIntersect(raytracer::Ray ray);
	};

}

#endif // !HEART_SHAPE_H