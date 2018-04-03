#ifndef RAY_H
#define RAY_H

#include "Vector.h"

namespace raytracer {

	class Ray {
		raytracer::Vector p0;
		raytracer::Vector v;

	public:
		Ray(raytracer::Vector p0 = raytracer::Vector(), raytracer::Vector v = raytracer::Vector());

		Ray& operator=(const Ray& other);

		Ray& P0(raytracer::Vector newP0);

		Ray& V(raytracer::Vector newV);

		Ray Reflect(raytracer::Vector point, raytracer::Vector normal);

		Ray Refract(raytracer::Vector point, raytracer::Vector normal, double n);

		raytracer::Vector GetP0();
		raytracer::Vector GetV();

	};

}

#endif // !RAY_H
