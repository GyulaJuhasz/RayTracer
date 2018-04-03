#include <cmath>

#include "MathUtils.h"
#include "Vector.h"

using raytracer::Vector;

Vector::Vector(double x, double y, double z, double h) : x(x), y(y), z(z), h(h) {}

double Vector::GetX() const { return x / h; }
double Vector::GetY() const { return y / h; }
double Vector::GetZ() const { return z / h; }
double Vector::GetH() const { return h; }

double Vector::GetRawX() const { return x; }
double Vector::GetRawY() const { return y; }
double Vector::GetRawZ() const { return z; }

Vector& Vector::X(double newX) {
	x = newX;
	return *this;
}
Vector& Vector::Y(double newY) {
	y = newY;
	return *this;
}
Vector& Vector::Z(double newZ) {
	z = newZ;
	return *this;
}
Vector& Vector::H(double newH) {
	h = newH;
	return *this;
}

void Vector::Set(double vX, double vY, double vZ, double vH) {
	x = vX;
	y = vY;
	z = vZ;
	h = vH;
}

double Vector::GetLenght() { return (double)sqrt((x*x + y*y + z*z) / (h*h)); }

Vector Vector::operator+(const Vector& v) { return Vector(GetX() + v.GetX(), GetY() + v.GetY(), GetZ() + v.GetZ(), 1.0); }
Vector Vector::operator-(const Vector& v) { return Vector(GetX() - v.GetX(), GetY() - v.GetY(), GetZ() - v.GetZ(), 1.0); }
Vector Vector::operator%(const Vector& v) { return Vector(GetY() * v.GetZ() - GetZ() * v.GetY(), GetZ() * v.GetX() - GetX() * v.GetZ(), GetX() * v.GetY() - GetY() * v.GetX()); }
Vector Vector::operator*(const double& c) { return Vector(c * GetX(), c * GetY(), c * GetZ(), 1.0); }
double Vector::operator*(const Vector& v) { return (GetX() * v.GetX() + GetY() * v.GetY() + GetZ() * v.GetZ()); }
Vector Vector::operator/(const double& c) { return Vector(GetX() / c, GetY() / c, GetZ() / c, 1.0); }

Vector& Vector::operator=(const Vector& v) {
	if (this != &v) {
		x = v.x;
		y = v.y;
		z = v.z;
		h = v.h;
	}
	return *this;
}

void Vector::Normalize() {
	double length = this->GetLenght();
	if (raytracer::math::IsZero(length)) {
		return;
	}

	double newX = (this->GetX()) / length;
	double newY = (this->GetY()) / length;
	double newZ = (this->GetZ()) / length;
	x = newX;
	y = newY;
	z = newZ;
	h = 1.0;
}

Vector Vector::GetNormalized() {
	Vector result(GetX(), GetY(), GetZ(), GetH());
	result.Normalize();
	return result;
}