#include "Ray.h"
#include "Config.h"
#include <cmath>

using raytracer::Ray;

Ray::Ray(raytracer::Vector p0, raytracer::Vector v) : p0(p0), v(v) {}

Ray& Ray::operator=(const Ray& other) {
	if (this != &other) {
		p0 = other.p0;
		v = other.v;
	}
	return *this;
}

Ray& Ray::P0(raytracer::Vector newP0) {
	p0 = newP0;
	return *this;
}

Ray& Ray::V(raytracer::Vector newV) {
	v = newV;
	return *this;
}

Ray Ray::Reflect(raytracer::Vector point, raytracer::Vector normal) {
	raytracer::Vector baseDir = v.GetNormalized();
	raytracer::Vector reflectedDir = baseDir - normal * (2.0 * (normal * baseDir));
	Ray reflected;
	reflected.P0(point + normal * EPSZILON_TOLAS).V(reflectedDir);
	return reflected;
}

Ray Ray::Refract(raytracer::Vector point, raytracer::Vector normal, double n) {
	raytracer::Vector baseDir = v.GetNormalized();
	double cosalfa = normal * (baseDir * (-1.0));
	double sinalfanegyzet = 1.0 - cosalfa*cosalfa;
	double cosbetanegyzet = 1.0 - sinalfanegyzet / (n*n);
	double cosbeta = sqrt(cosbetanegyzet);
	raytracer::Vector refractedDir = baseDir / n + normal * (cosalfa / n - cosbeta);
	Ray refracted;
	refracted.P0(point - normal * EPSZILON_TOLAS).V(refractedDir);
	return refracted;
}

raytracer::Vector Ray::GetP0() { return p0; }
raytracer::Vector Ray::GetV() { return v; }