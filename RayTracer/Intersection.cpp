#include "Intersection.h"

using raytracer::Intersection;

Intersection::Intersection():t(-1.0), index(-1) {}

Intersection& Intersection::Param(double param) {
	t = param;
	return *this;
}

double Intersection::GetParam() { return t; }

Intersection& Intersection::IntersectPoint(raytracer::Vector mp) {
	intersectPoint = mp;
	return *this;
}

raytracer::Vector Intersection::GetIntersectPoint() { return intersectPoint; }

Intersection& Intersection::Normal(raytracer::Vector norm) {
	normal = norm;
	return *this;
}

raytracer::Vector Intersection::GetNormal() { return normal; }

Intersection& Intersection::Index(int newIndex) {
	index = newIndex;
	return *this;
}

int Intersection::GetIndex() { return index; }

Intersection& Intersection::operator=(const Intersection& other) {
	if (this != &other) {
		intersectPoint = other.intersectPoint;
		normal = other.normal;
		t = other.t;
		index = other.index;
	}
	return *this;
}