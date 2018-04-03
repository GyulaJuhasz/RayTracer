#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H

#include "Config.h"
#include "Intersection.h"
#include "Object.h"
#include "Ray.h"
#include "Triangle.h"

namespace raytracer {

	class TriangleMesh : public raytracer::Object {
		raytracer::Triangle* triangles[MAX_HAROMSZOG_SZAM];
		int numOfTriangles;
	public:
		TriangleMesh();

		int GetNumberOfTriangles();

		void Clear();

		void AddTriangle(raytracer::Triangle *newTriangle);

		raytracer::Intersection Intersect(raytracer::Ray ray);

	};

}

#endif // !TRIANGLE_MESH_H
