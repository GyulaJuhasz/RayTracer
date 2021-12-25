#include <cmath>

#include "Ellipsoid.h"

using raytracer::Ellipsoid;

Ellipsoid::Ellipsoid() : raytracer::Object() {}

Ellipsoid& Ellipsoid::Center(raytracer::Vector newCenter) {
	center = newCenter;
	return *this;
}

raytracer::Vector Ellipsoid::GetCenter() { return center; }

Ellipsoid& Ellipsoid::A(double newA) {
	_A = newA;
	return *this;
}

double Ellipsoid::GetA() { return _A; }

Ellipsoid& Ellipsoid::B(double newB) {
	_B = newB;
	return *this;
}

double Ellipsoid::GetB() { return _B; }

Ellipsoid& Ellipsoid::C(double newC) {
	_C = newC;
	return *this;
}

double Ellipsoid::GetC() { return _C; }

raytracer::Intersection Ellipsoid::SpecificIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;

	result.Param(-1.0);

	Vector i(1.0, 0.0, 0.0);
	Vector j(0.0, 1.0, 0.0);
	Vector k(0.0, 0.0, 1.0);

	double px = ray.GetP0().GetX();
	double py = ray.GetP0().GetY();
	double pz = ray.GetP0().GetZ();
	double vx = ray.GetV().GetX();
	double vy = ray.GetV().GetY();
	double vz = ray.GetV().GetZ();
	double ex = center.GetX();
	double ey = center.GetY();
	double ez = center.GetZ();

	double eh1 = _B*_B * _C*_C;
	double eh2 = _A*_A * _C*_C;
	double eh3 = _A*_A * _B*_B;
	double R = _A*_A *_B*_B* _C*_C;

	double a = eh1*vx*vx + eh2*vy*vy + eh3*vz*vz;
	double b = 2.0 * (eh1 * (px*vx - vx*ex) + eh2 * (py*vy - vy*ey) + eh3 * (pz*vz - vz*ez));
	double c = eh1 * (px*px + ex*ex - 2.0 * px*ex) + eh2 * (py*py + ey*ey - 2.0 * py*ey) + eh3 * (pz*pz + ez*ez - 2.0 * pz*ez) - R;

	double discr = b*b - 4 * a*c;

	if (discr < 0.0) return result;

	double t1 = (-1.0*b - sqrt(discr)) / (2.0*a);
	double t2 = (-1.0*b + sqrt(discr)) / (2.0*a);

	if ((t1 > 0.0) || (t2 > 0.0)) {
		double t = (t1 > 0.0) ? t1 : t2;
		Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
		double ix = intersectPoint.GetX();
		double iy = intersectPoint.GetY();
		double iz = intersectPoint.GetZ();
		double nx = 2.0*(ix - ex) / (_A*_A);
		double ny = 2.0*(iy - ey) / (_B*_B);
		double nz = 2.0*(iz - ez) / (_C*_C);
		Vector normal = (i*nx + j*ny + k*nz).GetNormalized();
		result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
	}

	return result;
}