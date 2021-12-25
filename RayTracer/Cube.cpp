#include "Cube.h"

using raytracer::Cube;

Cube::Cube() {
	raytracer::Vector FLD(-2.0, -2.0, -2.0);
	raytracer::Vector FRD(2.0, -2.0, -2.0);
	raytracer::Vector FRU(2.0, -2.0, 2.0);
	raytracer::Vector FLU(-2.0, -2.0, 2.0);
	raytracer::Vector BLD(-2.0, 2.0, -2.0);
	raytracer::Vector BRD(2.0, 2.0, -2.0);
	raytracer::Vector BRU(2.0, 2.0, 2.0);
	raytracer::Vector BLU(-2.0, 2.0, 2.0);

	front.Corner1(FLD).Corner2(FRD).Corner3(FRU).Corner4(FLU);
	back.Corner1(BRD).Corner2(BLD).Corner3(BLU).Corner4(BRU);
	left.Corner1(BLD).Corner2(FLD).Corner3(FLU).Corner4(BLU);
	right.Corner1(FRD).Corner2(BRD).Corner3(BRU).Corner4(FRU);
	top.Corner1(FLU).Corner2(FRU).Corner3(BRU).Corner4(BLU);
	bottom.Corner1(FLD).Corner2(BLD).Corner3(BRD).Corner4(FRD);
}

raytracer::Intersection Cube::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	raytracer::Intersection actual;

	double tmin = 99999999999.0;
	actual = front.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	actual = back.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	actual = left.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	actual = right.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	actual = top.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	actual = bottom.Intersect(ray);
	if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
		tmin = actual.GetParam();
		result = actual;
	}
	return result;
}