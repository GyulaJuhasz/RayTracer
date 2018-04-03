#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Triangle : public raytracer::Object {
		raytracer::Vector vertex1;
		raytracer::Vector vertex2;
		raytracer::Vector vertex3;

	public:
		Triangle();

		Triangle& Vertex1(raytracer::Vector newVertex1);

		raytracer::Vector GetVertex1();

		Triangle& Vertex2(raytracer::Vector newVertex2);

		raytracer::Vector GetVertex2();

		Triangle& Vertex3(raytracer::Vector newVertex3);

		raytracer::Vector GetVertex3();

		raytracer::Intersection Intersect(raytracer::Ray ray);

		void Beallit(raytracer::Vector cs1, raytracer::Vector cs2, raytracer::Vector cs3);

	};

}

#endif // !TRIANGLE_H