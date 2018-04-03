#include "Triangle.h"

using raytracer::Triangle;

Triangle::Triangle() : raytracer::Object() {}

Triangle& Triangle::Vertex1(raytracer::Vector newVertex1) {
	vertex1 = newVertex1;
	return *this;
}

raytracer::Vector Triangle::GetVertex1() { return vertex1; }

Triangle& Triangle::Vertex2(raytracer::Vector newVertex2) {
	vertex2 = newVertex2;
	return *this;
}

raytracer::Vector Triangle::GetVertex2() { return vertex2; }

Triangle& Triangle::Vertex3(raytracer::Vector newVertex3) {
	vertex3 = newVertex3;
	return *this;
}

raytracer::Vector Triangle::GetVertex3() { return vertex3; }


raytracer::Intersection Triangle::Intersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	raytracer::Vector normal = ((vertex2 - vertex1) % (vertex3 - vertex1)).GetNormalized();

	raytracer::Vector p = ray.GetP0();
	raytracer::Vector v = ray.GetV();

	double t = (vertex1*normal - p*normal) / (v*normal);

	if (t > 0.0) {
		raytracer::Vector intersectPoint = p + v*t;

		double test1 = ((vertex2 - vertex1) % (intersectPoint - vertex1)) * normal;
		double test2 = ((vertex3 - vertex2) % (intersectPoint - vertex2)) * normal;
		double test3 = ((vertex1 - vertex3) % (intersectPoint - vertex3)) * normal;

		bool good1 = (test1 > 0.0);
		bool good2 = (test2 > 0.0);
		bool good3 = (test3 > 0.0);

		if (good1 && good2 && good3) {
			result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
		}
	}

	return result;
}

void Triangle::Beallit(raytracer::Vector cs1, raytracer::Vector cs2, raytracer::Vector cs3) {
	vertex1 = cs1;
	vertex2 = cs2;
	vertex3 = cs3;
}
