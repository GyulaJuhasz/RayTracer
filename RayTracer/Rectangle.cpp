#include "Rectangle.h"
#include <cmath>

using raytracer::Rectangle;

Rectangle::Rectangle() : raytracer::Object() {}

Rectangle& Rectangle::Corner1(raytracer::Vector newCorner1) {
	corner1 = newCorner1;
	return *this;
}

raytracer::Vector Rectangle::GetCorner1() { return corner1; }

Rectangle& Rectangle::Corner2(raytracer::Vector newCorner2) {
	corner2 = newCorner2;
	return *this;
}

raytracer::Vector Rectangle::GetCorner2() { return corner2; }

Rectangle& Rectangle::Corner3(raytracer::Vector newCorner3) {
	corner3 = newCorner3;
	return *this;
}

raytracer::Vector Rectangle::GetCorner3() { return corner3; }

Rectangle& Rectangle::Corner4(raytracer::Vector newCorner4) {
	corner4 = newCorner4;
	return *this;
}

raytracer::Vector Rectangle::GetCorner4() { return corner4; }

double Rectangle::GetSizeX() { return (corner2 - corner1).GetLenght(); }
double Rectangle::GetSizeY() { return (corner3 - corner1).GetLenght(); }

raytracer::Intersection Rectangle::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	raytracer::Vector normal = ((corner2 - corner1) % (corner3 - corner1)).GetNormalized();

	raytracer::Vector p = ray.GetP0();
	raytracer::Vector v = ray.GetV();

	double t = (corner1*normal - p*normal) / (v*normal);

	if (t >= 0.0) {
		raytracer::Vector intersectPoint = p + v*t;

		double test1 = ((corner2 - corner1) % (intersectPoint - corner1)) * normal;
		double test2 = ((corner3 - corner2) % (intersectPoint - corner2)) * normal;
		double test3 = ((corner4 - corner3) % (intersectPoint - corner3)) * normal;
		double test4 = ((corner1 - corner4) % (intersectPoint - corner4)) * normal;

		bool good1 = (test1 > 0.0);
		bool good2 = (test2 > 0.0);
		bool good3 = (test3 > 0.0);
		bool good4 = (test4 > 0.0);

		if (good1 && good2 && good3 && good4) {
			result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
		}
	}

	return result;
}

void Rectangle::AddPhoton(raytracer::Color intensity, raytracer::Vector position) {
	if (photonMap != NULL) {
		raytracer::Vector a = (corner2 - corner1).GetNormalized();
		raytracer::Vector c = position - corner1;
		double e = c.GetLenght();

		double photonCellSizeX = GetSizeX() / FOTONTERKEP_X;
		double photonCellSizeY = GetSizeY() / FOTONTERKEP_Y;

		int photonX = 0;
		int	photonY = 0;

		if (e > 0.0) {
			c.Normalize();
			double cosalfa = c*a;
			double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

			double pointX = e*cosalfa;
			double pointY = e*sinalfa;

			while (pointX >= photonCellSizeX) {
				photonX++;
				pointX -= photonCellSizeX;
			}

			while (pointY >= photonCellSizeY) {
				photonY++;
				pointY -= photonCellSizeY;
			}
		}

		photonMap->AddPhoton(photonX, photonY, intensity);
	}

}

raytracer::Color Rectangle::GetPhoton(raytracer::Vector position) {
	raytracer::Color photon = raytracer::Object::GetPhoton(position);

	if (photonMap != NULL) {
		raytracer::Vector a = (corner2 - corner1).GetNormalized();
		raytracer::Vector c = position - corner1;
		double e = c.GetLenght();

		double photonCellSizeX = GetSizeX() / FOTONTERKEP_X;
		double photonCellSizeY = GetSizeY() / FOTONTERKEP_Y;

		int photonX = 0;
		int	photonY = 0;

		if (e > 0.0) {
			c.Normalize();
			double cosalfa = c*a;
			double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

			double pointX = e*cosalfa;
			double pointY = e*sinalfa;

			while (pointX >= photonCellSizeX) {
				photonX++;
				pointX -= photonCellSizeX;
			}

			while (pointY >= photonCellSizeY) {
				photonY++;
				pointY -= photonCellSizeY;
			}
		}

		photon = photonMap->GetPhoton(photonX, photonY);
	}

	return photon;
}

raytracer::Color Rectangle::GetProceduralColor(raytracer::Vector position) {
	raytracer::Color result = raytracer::Object::GetProceduralColor(position);

	if (TEXTURA) {
		int i = 0;
		int j = 0;

		raytracer::Vector a = (corner2 - corner1).GetNormalized();

		double unitX = GetSizeX() / TEXTURA_X;
		double unitY = GetSizeY() / TEXTURA_Y;

		raytracer::Vector c = position - corner1;
		double e = c.GetLenght();

		if (e > 0.0) {
			c.Normalize();
			double cosalfa = c*a;
			double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

			double diffX = e*cosalfa;
			double diffY = e*sinalfa;

			while (diffX > unitX) {
				diffX -= unitX;
				i++;
			}
			while (diffY > unitY) {
				diffY -= unitY;
				j++;
			}
		}

		double factor = sqrt((double)(i^j)) + sqrt((double)(j^i));
		int multiplier = (int)(floor(factor)) % 5;

		result = result * (0.1 + 0.1*multiplier);
	}

	return result;
}
