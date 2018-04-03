#include "HeartShape.h"
#include "MathUtils.h"
#include "Polinom.h"
#include "RealRoots.h"
#include "Vector.h"

using raytracer::HeartShape;

HeartShape::HeartShape() {
	minX = 99999999.0;
	maxX = -99999999.0;
	minY = 99999999.0;
	maxY = -99999999.0;
	minZ = 99999999.0;
	maxZ = -99999999.0;
}

raytracer::Intersection HeartShape::Intersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	raytracer::Intersection test;
	test = cube.Intersect(ray);
	double tmin = test.GetParam();

	if (tmin < 0.0) return result;

	raytracer::Vector newP0 = test.GetIntersectPoint() + ray.GetV() * 0.001;
	raytracer::Ray ray2(newP0, ray.GetV());

	raytracer::Intersection test2 = cube.Intersect(ray2);
	double ttemp = test2.GetParam();

	if (ttemp < 0.0) return result;

	raytracer::Vector testIntersect = test2.GetIntersectPoint();

	double tmax = -1.0;

	double tix = testIntersect.GetX();
	double tiy = testIntersect.GetY();
	double tiz = testIntersect.GetZ();

	double px = ray.GetP0().GetX();
	double py = ray.GetP0().GetY();
	double pz = ray.GetP0().GetZ();

	double vx = ray.GetV().GetX();
	double vy = ray.GetV().GetY();
	double vz = ray.GetV().GetZ();

	if (!raytracer::math::IsZero(vx)) {
		tmax = (tix - px) / vx;
	}
	else if (!raytracer::math::IsZero(vy)) {
		tmax = (tiy - py) / vy;
	}
	else if (!raytracer::math::IsZero(vy)) {
		tmax = (tiz - pz) / vz;
	}

	if (tmax < 0.0) return result;

	//double guess = test.GetParam();

	double c1 = vx * vx + 2.25 * vy * vy + vz * vz;
	double c2 = 2 * px * vx + 4.5 * py * vy + 2 * pz * vz;
	double c3 = px * px + 2.25 * py * py + pz * pz - 1;

	double c4 = c1 * c1;
	double c5 = 2 * c1 * c2;
	double c6 = c2 * c2 + 2 * c1 * c3;
	double c7 = 2 * c2 * c3;
	double c8 = c3 * c3;

	double c9 = c1 * c4;
	double c10 = c1 * c5 + c2 * c4;
	double c11 = c1 * c6 + c2 * c5 + c3 * c4;
	double c12 = c1 * c7 + c2 * c6 + c3 * c5;
	double c13 = c1 * c8 + c2 * c7 + c3 * c6;
	double c14 = c2 * c8 + c3 * c7;
	double c15 = c3 * c8;

	double c16 = vy * vy;
	double c17 = 2 * py * vy;
	double c18 = py * py;

	double c19 = vz * vz * vz;
	double c20 = 3 * pz * vz * vz;
	double c21 = 3 * pz * pz * vz;
	double c22 = pz * pz * pz;

	double c23 = c16 * c19;
	double c24 = c16 * c20 + c17 * c19;
	double c25 = c16 * c21 + c17 * c20 + c18 * c19;
	double c26 = c16 * c22 + c17 * c21 + c18 * c20;
	double c27 = c17 * c22 + c18 * c21;
	double c28 = c18 * c22;

	double c29 = vx * vx;
	double c30 = 2 * px * vx;
	double c31 = px * px;

	double c32 = c29 * c19;
	double c33 = c29 * c20 + c30 * c19;
	double c34 = c29 * c21 + c30 * c20 + c31 * c19;
	double c35 = c29 * c22 + c30 * c21 + c31 * c20;
	double c36 = c30 * c22 + c31 * c21;
	double c37 = c31 * c22;

	double constant1 = -9.0 / 80;

	double constant2 = -1.0;

	double a6 = c9;
	double a5 = c10 + constant1 * c23 + constant2 * c32;
	double a4 = c11 + constant1 * c24 + constant2 * c33;
	double a3 = c12 + constant1 * c25 + constant2 * c34;
	double a2 = c13 + constant1 * c26 + constant2 * c35;
	double a1 = c14 + constant1 * c27 + constant2 * c36;
	double a0 = c15 + constant1 * c28 + constant2 * c37;

	raytracer::Polinom heartPoli(6);

	heartPoli.SetCoefficient(6, a6);
	heartPoli.SetCoefficient(5, a5);
	heartPoli.SetCoefficient(4, a4);
	heartPoli.SetCoefficient(3, a3);
	heartPoli.SetCoefficient(2, a2);
	heartPoli.SetCoefficient(1, a1);
	heartPoli.SetCoefficient(0, a0);

	//RealRoots roots = heartPoli.getRootsWithNewton(10.0);
	raytracer::RealRoots roots = heartPoli.getRootsWithBisection(heartPoli, tmin, tmax, NULL);

	/*
	ofstream newtonFile("newton.txt", ios::out | ios::app);
	if ( newtonFile.is_open() ) {
	newtonFile << "---------------------------------------" << endl;
	newtonFile << "---------------------------------------" << endl;
	newtonFile << heartPoli << " solutions: " << endl;
	newtonFile << "---------------------------------------" << endl;
	newtonFile << roots << endl;
	newtonFile << "---------------------------------------" << endl;
	newtonFile << "---------------------------------------" << endl;
	newtonFile.close();
	}
	*/

	int number = roots.GetNumberOfRoots();

	if (number == 0) {
		return result;
	}

	double t = roots.GetMinPositiveRootIfAny();

	if (t <= 0.0) {
		return result;
	}

	raytracer::Vector intersectPoint = ray.GetP0() + ray.GetV() * t;

	double ix = intersectPoint.GetX();
	double iy = intersectPoint.GetY();
	double iz = intersectPoint.GetZ();

	if (ix < minX) minX = ix;
	if (ix > maxX) maxX = ix;
	if (iy < minY) minY = iy;
	if (iy > maxY) maxY = iy;
	if (iz < minZ) minZ = iz;
	if (iz > maxZ) maxZ = iz;

	double n1 = ix * ix + 2.25f * iy * iy + iz * iz - 1.0;
	double n2 = n1 * n1;

	double parcX = 6 * ix * n2 - 2 * ix * iz * iz * iz;
	double parcY = 13.5 * iy * n2 - (9.0 / 40) * iy * iz * iz * iz;
	double parcZ = 6 * iz * n2 - 3 * ix * ix * iz * iz - (27.0 / 80) * iy * iy * iz * iz;

	/*
	double dxdy = parcX / parcY;
	double dxdz = parcX / parcZ;
	double dydx = parcY / parcX;
	double dydz = parcY / parcZ;
	double dzdx = parcZ / parcX;
	double dzdy = parcZ / parcY;

	double nx = 2 * n2 * (2 * ix + 2.25 * 2 * iy * dydx + 2 * iz * dzdx ) - (2 * ix * iz * iz * iz + ix * ix * 3 * iz * iz * dzdx ) - (9.0/80) * (2*iy*dydx * iz*iz*iz + iy*iy*3*iz*iz*dzdx );
	double ny = 2 * n2 * (2 * ix * dxdy + 2.25 * 2 * iy + 2 * iz * dzdy) - (2 * ix * dxdy * iz * iz * iz + ix * ix * 3 * iz * iz * dzdy ) - ( 9.0/80 ) * (2* iy * iz * iz * iz + iy * iy * 3 * iz * iz * dzdy );
	double nz = 2 * n2 * (2 * ix * dxdz + 2.25 * 2 * iy * dydz + 2 * iz ) - (2 * ix * dxdz * iz * iz * iz + ix * ix * 3 * iz * iz ) - (9.0 / 80 ) * (2 * iy * dydz * iz * iz * iz + iy * iy * 3 * iz * iz);

	Vector normal(nx, ny, nz);
	*/
	raytracer::Vector normal(parcX, parcY, parcZ);
	normal.Normalize();
	result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
	return result;
}
