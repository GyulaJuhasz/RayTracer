#include <cstddef>

#include "CompositeObject.h"

using raytracer::CompositeObject;

CompositeObject::CompositeObject() : raytracer::Object() {
	RemoveAllParts();
}

void CompositeObject::RemoveAllParts() {
	numberOfParts = 0;
	for (int i = 0; i < MAX_NUMBER_OF_COMPOSITE_OBJECT_PARTS; i++) {
		parts[i] = NULL;
	}
}

CompositeObject& CompositeObject::AddPart(raytracer::Object* newPart) {
	if (numberOfParts < MAX_NUMBER_OF_COMPOSITE_OBJECT_PARTS && newPart != NULL) {
		parts[numberOfParts++] = newPart;
	}
	return *this;
}

raytracer::Intersection CompositeObject::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	raytracer::Intersection actualIntersection;

	double closestIntersectionParam = 99999999999.0;
	if (numberOfParts > 0) {
		for (int i = 0; i < numberOfParts; i++) {
			actualIntersection = parts[i]->Intersect(ray);
			double intersectionParam = actualIntersection.GetParam();
			if ((intersectionParam > 0.0) && (intersectionParam < closestIntersectionParam)) {
				closestIntersectionParam = intersectionParam;
				result = actualIntersection;
				result.Index(i);
			}
		}
	}
	return result;

}