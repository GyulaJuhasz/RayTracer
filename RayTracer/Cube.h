#ifndef CUBE_H

#define CUBE_H

#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Rectangle.h"

namespace raytracer {

	class Cube : public raytracer::Object {
		raytracer::Rectangle front;
		raytracer::Rectangle back;
		raytracer::Rectangle left;
		raytracer::Rectangle right;
		raytracer::Rectangle top;
		raytracer::Rectangle bottom;

	public:
		Cube();

		raytracer::Intersection Intersect(raytracer::Ray ray);

	};

}

#endif // !CUBE_H