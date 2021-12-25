#include "Cylinder.h"
#include "MathUtils.h"

using raytracer::Cylinder;

Cylinder::Cylinder() : raytracer::Object() { solid = false; }

Cylinder& Cylinder::BasePoint(raytracer::Vector newBasePoint) {
	basePoint = newBasePoint;
	return *this;
}

raytracer::Vector Cylinder::GetBasePoint() { return basePoint; }

Cylinder& Cylinder::Radius(double newRadius) {
	radius = newRadius;
	return *this;
}

double Cylinder::GetRadius() { return radius; }

Cylinder& Cylinder::Height(double newHeight) {
	height = newHeight;
	return *this;
}

double Cylinder::GetHeight() { return height; }

Cylinder& Cylinder::Solid(bool newSolid) {
	solid = newSolid;
	return *this;
}

bool Cylinder::IsSolid() { return solid; }

raytracer::Intersection Cylinder::Intersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	result.Param(-1.0);

	double px = ray.GetP0().GetX();
	double py = ray.GetP0().GetY();
	double vx = ray.GetV().GetX();
	double vy = ray.GetV().GetY();
	double hx = basePoint.GetX();
	double hy = basePoint.GetY();
	double r = radius;

	double a = vx*vx + vy*vy;
	double b = 2.0*(px*vx + py*vy - vx*hx - vy*hy);
	double c = px*px + py*py + hx*hx + hy*hy - 2.0*(px*hx + py*hy) - r*r;
	double discr = b*b - 4 * a*c;
	if (discr < 0.0) return result;

	double bottomHeight = basePoint.GetZ();
	double topHeight = bottomHeight + height;
	if (raytracer::math::IsZero(discr)) {
		double t = (-1.0*b) / (2.0*a);
		Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
		Vector helper;
		double z = intersectPoint.GetZ();
		if ((z < bottomHeight) ||
			(z > topHeight)
			|| (t < 0.0)) {
			return result;
		}
		else {
			helper = intersectPoint;
			helper.Z(bottomHeight);
			Vector normal = (helper - basePoint).GetNormalized();
			result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
			return result;
		}
	}
	else {
		double t1 = (-1.0*b - sqrt(discr)) / (2.0*a);
		double t2 = (-1.0*b + sqrt(discr)) / (2.0*a);

		Vector intersectPoint1 = ray.GetP0() + ray.GetV()*t1;
		Vector intersectPoint2 = ray.GetP0() + ray.GetV()*t2;
		Vector helper;
		double z1 = intersectPoint1.GetZ();
		double z2 = intersectPoint2.GetZ();

		if ((z1 >= bottomHeight) && (z1 <= topHeight) && (t1 > 0.0)) {
			// Elsõ metszéspont magassága stimmel
			helper = intersectPoint1;
			helper.Z(bottomHeight);
			Vector normal = (helper - basePoint).GetNormalized();
			result.Param(t1).IntersectPoint(intersectPoint1).Normal(normal);
			return result;
		}

		if (solid) {
			// Tömör hengernél ellenõrizni kell a tetõ és alaplap metszését
			double pz = ray.GetP0().GetZ();
			double vz = ray.GetV().GetZ();

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

			if (firstHeight == ABOVE && (secondHeight == CORRECT || secondHeight == BELOW)) {
				// Tetõlap metszve

				// pz + vz*t = teto => t = (teto - pz) / vz
				double t = (topHeight - pz) / vz;
				Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
				Vector normal(0.0, 0.0, 1.0);
				result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
				return result;
			}

			if (firstHeight == BELOW && (secondHeight == CORRECT || secondHeight == ABOVE)) {
				// Alaplap metszve

				// pz + vz*t = alapraytracer::Colort => t = (alapraytracer::Colort - pz) / vz
				double t = (bottomHeight - pz) / vz;
				Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
				Vector normal(0.0, 0.0, -1.0);
				result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
				return result;
			}

		}

		// Nem tömör hengernél eddigi módszer
		if ((z2 > bottomHeight) && (z2 < topHeight) && (t2 > 0.0)) {
			helper = intersectPoint2;
			helper.Z(bottomHeight);
			Vector normal = (helper - basePoint).GetNormalized();
			normal = normal * -1.0;
			result.Param(t2).IntersectPoint(intersectPoint2).Normal(normal);
			return result;
		}

		return result;
	}
}
