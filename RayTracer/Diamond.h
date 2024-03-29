#ifndef DIAMOND_H
#define DIAMOND_H

#include "Ellipsoid.h"
#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "TriangleMesh.h"

namespace raytracer {

	class Diamond : public raytracer::Object {
		raytracer::Ellipsoid *boundingVolume;
		raytracer::TriangleMesh *diamondBody;
	public:

		Diamond();

		Diamond& DiamondBody(raytracer::TriangleMesh *newBody);

		raytracer::TriangleMesh* GetDiamondBody();

	protected:
		raytracer::Intersection SpecificIntersect(raytracer::Ray ray);
	};

}

#endif // !DIAMOND_H