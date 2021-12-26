#ifndef COMPOSITE_OBJECT_H
#define COMPOSITE_OBJECT_H

#include "Config.h"
#include "Intersection.h"
#include "Object.h"
#include "Ray.h"

namespace raytracer {

	class CompositeObject : public raytracer::Object {
		raytracer::Object* parts[MAX_NUMBER_OF_COMPOSITE_OBJECT_PARTS];
		int numberOfParts;

	public:
		CompositeObject();

		void RemoveAllParts();

		CompositeObject& AddPart(raytracer::Object* newPart);

	protected:
		Intersection SpecificIntersect(raytracer::Ray ray);
	};
}

#endif // !COMPOSITE_OBJECT_H