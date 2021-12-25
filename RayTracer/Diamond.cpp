#include <cstddef>

#include "Diamond.h"

using raytracer::Diamond;

Diamond::Diamond() {
	boundingVolume = NULL;
	diamondBody = NULL;
}

Diamond& Diamond::DiamondBody(raytracer::TriangleMesh* newBody) {
	diamondBody = newBody;
	return *this;
}

raytracer::TriangleMesh* Diamond::GetDiamondBody() { return diamondBody; }

raytracer::Intersection Diamond::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	if (diamondBody != NULL) {
		result = diamondBody->Intersect(ray);
	}

	return result;
}