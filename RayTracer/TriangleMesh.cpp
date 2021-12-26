#include <cstddef>

#include "TriangleMesh.h"

using raytracer::TriangleMesh;

TriangleMesh::TriangleMesh() : raytracer::Object() {
	Clear();
}

int TriangleMesh::GetNumberOfTriangles() { return numOfTriangles; }

void TriangleMesh::Clear() {
	numOfTriangles = 0;
	for (int i = 0; i < MAX_HAROMSZOG_SZAM; i++) triangles[i] = NULL;
}

void TriangleMesh::AddTriangle(raytracer::Triangle *newTriangle) {
	if (numOfTriangles < MAX_HAROMSZOG_SZAM && newTriangle != NULL) {
		triangles[numOfTriangles++] = newTriangle;
	}
}

raytracer::Intersection TriangleMesh::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	raytracer::Intersection actual;

	double tmin = 99999999999.0;
	if (numOfTriangles > 0) {
		for (int i = 0; i < numOfTriangles; i++) {
			actual = triangles[i]->Intersect(ray);
			if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
				tmin = actual.GetParam();
				result = actual;
				result.Index(i);
			}
		}
	}
	return result;
}