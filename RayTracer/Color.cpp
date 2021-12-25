#include "Color.h"

using raytracer::Color;

Color::Color(double r, double g, double b) : r(r), g(g), b(b) {}

double Color::GetR() const { return r; }
double Color::GetG() const { return g; }
double Color::GetB() const { return b; }

Color& Color::R(double newR) {
	r = newR;
	return *this;
}

Color& Color::G(double newG) {
	g = newG;
	return *this;
}

Color& Color::B(double newB) {
	b = newB;
	return *this;
}

void Color::Set(double cR, double cG, double cB) {
	r = cR;
	g = cG;
	b = cB;
}

Color& Color::operator=(const Color& other) {
	if (this != &other) {
		r = other.r;
		g = other.g;
		b = other.b;
	}
	return *this;
}

Color Color::operator+(const Color& other) { return Color(GetR() + other.GetR(), GetG() + other.GetG(), GetB() + other.GetB()); }
Color Color::operator-(const Color& other) { return Color(GetR() - other.GetR(), GetG() - other.GetG(), GetB() - other.GetB()); }
Color Color::operator*(const double& c) { return Color(GetR() * c, GetG() * c, GetB() * c); }
Color Color::operator/(const double& c) { return Color(GetR() / c, GetG() / c, GetB() / c); }
Color Color::operator*(const Color& other) { return Color(GetR() * other.GetR(), GetG() * other.GetG(), GetB() * other.GetB()); }
Color Color::operator/(const Color& other) { return Color(GetR() / other.GetR(), GetG() / other.GetG(), GetB() / other.GetB()); }