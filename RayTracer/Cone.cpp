#include "Cone.h"
#include "MathUtils.h"

using raytracer::Cone;

Cone::Cone() : raytracer::Object(), baseCenter(Vector(0.0, 0.0, 0.0)), baseRadius(1.0), height(1.0), isSolid(true) {
	CalculateValues();
}

Cone& Cone::BaseCenter(raytracer::Vector newBaseCenter) {
	baseCenter = newBaseCenter;
	CalculateValues();
	return *this;
}

Cone& Cone::BaseRadius(double newBaseRadius) {
	baseRadius = newBaseRadius;
	CalculateValues();
	return *this;
}

Cone& Cone::Height(double newHeight) {
	height = newHeight;
	CalculateValues();
	return *this;
}

Cone& Cone::IsSolid(bool newIsSolid) {
	isSolid = newIsSolid;
	return *this;
}

void Cone::CalculateValues() {
	vertex = baseCenter + Vector(0, 0, height);
	baseRadiusPerHeight = baseRadius / height;
	double denominator = sqrt(baseRadius * baseRadius + height * height);
	sinAlpha = baseRadius / denominator;
	cosAlpha = height / denominator;
	oneMinusCosAlpha = 1.0 - cosAlpha;
}

raytracer::Intersection Cone::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	result.Param(-1.0);

	raytracer::Vector p0 = ray.GetP0();
	raytracer::Vector v = ray.GetV();

	double p0x = p0.GetX();
	double p0y = p0.GetY();
	double p0z = p0.GetZ();
	double vx = v.GetX();
	double vy = v.GetY();
	double vz = v.GetZ();

	double cx = vertex.GetX();
	double cy = vertex.GetY(); 
	double cz = vertex.GetZ();

	double a = vx * vx + vy * vy - baseRadiusPerHeight * baseRadiusPerHeight * vz * vz;
	double b = 2.0 * (p0x * vx - vx * cx + p0y * vy - vy * cy - baseRadiusPerHeight * baseRadiusPerHeight * (p0z * vz - cz * vz));
	double c = p0x * p0x - 2 * p0x * cx + cx * cx + p0y * p0y - 2 * p0y * cy + cy * cy - baseRadiusPerHeight * baseRadiusPerHeight * (cz * cz + p0z * p0z - 2 * cz * p0z);
	double discriminant = b * b - 4 * a * c;
	if (discriminant < 0.0) {
		// No root, no intersections
		return result;
	}

	double bottomHeight = baseCenter.GetZ();
	double topHeight = cz;
	if (raytracer::math::IsZero(discriminant)) {
		// One root, one intersection
		double t = (-1.0 * b) / (2.0 * a);
		Vector intersectPoint = p0 + v * t;
		double z = intersectPoint.GetZ();
		if ((z < bottomHeight) || (z > topHeight) || (t < 0.0)) {
			// Intersection is not within height boundaries
			return result;
		}
		result.Param(t).IntersectPoint(intersectPoint).Normal(GetSurfaceNormal(intersectPoint));
		return result;
	}

	// Two roots, two possible intersections
	double squareRootOfDiscriminant = sqrt(discriminant);
	double t1 = (-1.0 * b - squareRootOfDiscriminant) / (2.0 * a);
	double t2 = (-1.0 * b + squareRootOfDiscriminant) / (2.0 * a);

	Vector intersectPoint1 = p0 + v * t1;
	Vector intersectPoint2 = p0 + v * t2;
	double z1 = intersectPoint1.GetZ();
	double z2 = intersectPoint2.GetZ();

	if ((z1 >= bottomHeight) && (z1 <= topHeight) && (t1 > 0.0)) {
		// First intersection is within height boundaries
		result.Param(t1).IntersectPoint(intersectPoint1).Normal(GetSurfaceNormal(intersectPoint1));
		return result;
	}

	if (isSolid) {
		// Solid cone has a base that could be intersected
		int ABOVE = 1;
		int CORRECT = 2;
		int BELOW = 3;

		int firstHeight, secondHeight;

		if (z1 > topHeight) {
			firstHeight = ABOVE;
		}
		else if (z1 < bottomHeight) {
			firstHeight = BELOW;
		}
		else {
			firstHeight = CORRECT;
		}

		if (z2 > topHeight) {
			secondHeight = ABOVE;
		}
		else if (z2 < bottomHeight) {
			secondHeight = BELOW;
		}
		else {
			secondHeight = CORRECT;
		}

		if (firstHeight == BELOW && (secondHeight == CORRECT || secondHeight == ABOVE)) {
			// Base intersected

			// p0z + vz*t = bottomHeight => t = (bottomHeight - p0z) / vz
			double t = (bottomHeight - p0z) / vz;
			Vector intersectPoint = p0 + v * t;
			Vector normal(0.0, 0.0, -1.0);
			result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
			return result;
		}
	}

	if ((z2 > bottomHeight) && (z2 < topHeight) && (t2 > 0.0)) {
		// Second intersection is reaching the surface from the inside
		raytracer::Vector normal = GetSurfaceNormal(intersectPoint2) * -1.0;
		result.Param(t2).IntersectPoint(intersectPoint2).Normal(normal);
		return result;
	}

	return result;
}

raytracer::Vector Cone::GetSurfaceNormal(raytracer::Vector surfacePoint) {
	raytracer::Vector centerPoint = raytracer::Vector(baseCenter.GetX(), baseCenter.GetY(), surfacePoint.GetZ());
	raytracer::Vector r = surfacePoint - centerPoint;
	raytracer::Vector upVector = vertex - centerPoint;
	raytracer::Vector rotationAxis = (r % upVector).GetNormalized();
	raytracer::Vector normal = r * cosAlpha + (rotationAxis % r) * sinAlpha + rotationAxis * (rotationAxis * r) * oneMinusCosAlpha;
	return normal.GetNormalized();
}