#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "Color.h"
#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Rectangle : public raytracer::Object {
		raytracer::Vector corner1;
		raytracer::Vector corner2;
		raytracer::Vector corner3;
		raytracer::Vector corner4;
	public:
		Rectangle();

		Rectangle& Corner1(raytracer::Vector newCorner1);

		raytracer::Vector GetCorner1();

		Rectangle& Corner2(raytracer::Vector newCorner2);

		raytracer::Vector GetCorner2();

		Rectangle& Corner3(raytracer::Vector newCorner3);

		raytracer::Vector GetCorner3();

		Rectangle& Corner4(raytracer::Vector newCorner4);

		raytracer::Vector GetCorner4();

		double GetSizeX();
		double GetSizeY();

		Intersection Intersect(raytracer::Ray ray);

		virtual void AddPhoton(raytracer::Color intensity, raytracer::Vector position);

		virtual raytracer::Color GetPhoton(raytracer::Vector position);

		virtual raytracer::Color GetProceduralColor(raytracer::Vector position);

	};

}

#endif // !RECTANGLE_H
