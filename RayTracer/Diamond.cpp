#include <cstddef>

#include "Diamond.h"

using raytracer::Diamond;

Diamond::Diamond() {
	boundingVolume = NULL;
	diamondBody = NULL;
}

Diamond& Diamond::BoundingVolume(raytracer::Ellipsoid *newBounding) {
	boundingVolume = newBounding;
	return *this;
}

raytracer::Ellipsoid* Diamond::GetBoundingVolume() { return boundingVolume; }

Diamond& Diamond::DiamondBody(raytracer::TriangleMesh *newBody) {
	diamondBody = newBody;
	return *this;
}

raytracer::TriangleMesh* Diamond::GetDiamondBody() { return diamondBody; }

raytracer::Intersection Diamond::Intersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	if (boundingVolume != NULL && diamondBody != NULL) {
		raytracer::Intersection boundingIntersection = boundingVolume->Intersect(ray);

		if (boundingIntersection.GetParam() >= 0.0) {
			result = diamondBody->Intersect(ray);
		}
	}

	return result;
}