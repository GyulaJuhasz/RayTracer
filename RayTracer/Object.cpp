#include "Object.h"
#include <cmath>

using raytracer::Object;

Object::Object() {
	isDiffuse = false;
	isReflective = false;
	isRefractive = false;
	kd.R(0.0).G(0.0).B(0.0);
	ks.R(0.0).G(0.0).B(0.0);
	ka.R(0.0).G(0.0).B(0.0);
	shine = 0.0;
	F0.R(0.0).G(0.0).B(0.0);
	n.R(0.0).G(0.0).B(0.0);
	k.R(0.0).G(0.0).B(0.0);
	proceduralMode = ProceduralMode::NONE;
	photonMap = NULL;
}

Object& Object::Diffuse(bool diffuse) {
	isDiffuse = diffuse;
	return *this;
}

bool Object::IsDiffuse() { return isDiffuse; }

Object& Object::Reflective(bool reflective) {
	isReflective = reflective;
	return *this;
}

bool Object::IsReflective() { return isReflective; }

Object& Object::Refractive(bool refractive) {
	isRefractive = refractive;
	return *this;
}

bool Object::IsRefractive() { return isRefractive; }

Object& Object::Kd(raytracer::Color newKd) {
	kd = newKd;
	return *this;
}

raytracer::Color Object::GetKd() { return kd; }

Object& Object::Ks(raytracer::Color newKs) {
	ks = newKs;
	return *this;
}

raytracer::Color Object::GetKs() { return ks; }

Object& Object::Ka(raytracer::Color newKa) {
	ka = newKa;
	return *this;
}

raytracer::Color Object::GetKa() { return ka; }

Object& Object::Shine(double newShine) {
	shine = newShine;
	return *this;
}

double Object::GetShine() { return shine; }

void Object::calculateFresnel() {
	F0
		.R(((n.GetR() - 1.0)*(n.GetR() - 1.0) + k.GetR()*k.GetR()) / ((n.GetR() + 1.0)*(n.GetR() + 1.0) + k.GetR()*k.GetR()))
		.G(((n.GetG() - 1.0)*(n.GetG() - 1.0) + k.GetG()*k.GetG()) / ((n.GetG() + 1.0)*(n.GetG() + 1.0) + k.GetG()*k.GetG()))
		.B(((n.GetB() - 1.0)*(n.GetB() - 1.0) + k.GetB()*k.GetB()) / ((n.GetB() + 1.0)*(n.GetB() + 1.0) + k.GetB()*k.GetB()));
}

Object& Object::N(raytracer::Color newN) {
	n = newN;
	calculateFresnel();
	return *this;
}

raytracer::Color Object::GetN() { return n; }

Object& Object::K(raytracer::Color newK) {
	k = newK;
	calculateFresnel();
	return *this;
}

Object& Object::BoundingObject(raytracer::Object *newBoundingObject) {
	boundingObject = newBoundingObject;
	return *this;
}

raytracer::Color Object::GetK() { return k; }

raytracer::Color Object::GetFresnel(raytracer::Vector N, raytracer::Vector V) {
	double cosalfa = fabs(N*V);
	raytracer::Color result;
	raytracer::Color a;
	raytracer::Color b;
	double factor = pow((1.0 - cosalfa), 5.0);
	a = F0;
	b.R(1.0 - F0.GetR()).G(1.0 - F0.GetG()).B(1.0 - F0.GetB());
	result = a + b*factor;
	return result;
}

Object& Object::ProcMode(ProceduralMode newMode) {
	proceduralMode = newMode;
	return *this;
}

Object::ProceduralMode Object::GetProcMode() { return proceduralMode; }

raytracer::Intersection Object::Intersect(raytracer::Ray ray) {
	if (boundingObject != NULL) {
		raytracer::Intersection boundingObjectIntersection = boundingObject->Intersect(ray);
		if (boundingObjectIntersection.GetParam() < 0.0) {
			return boundingObjectIntersection;
		}
	}
	return SpecificIntersect(ray);
}

void Object::SetPhotonMap(raytracer::PhotonMap *newMap) {
	photonMap = newMap;
}

raytracer::PhotonMap* Object::GetPhotonMap() {
	return photonMap;
}

void Object::WritePhotonMapToFile(const char* fileName) {
	if (photonMap != NULL) photonMap->WriteToFile(fileName);
}

void Object::ReadPhotonMapFromFile(const char* fileName) {
	if (photonMap != NULL) photonMap->ReadFromFile(fileName);
}

void Object::AddPhoton(raytracer::Color intensity, raytracer::Vector position) {};

raytracer::Color Object::GetPhoton(raytracer::Vector position) {
	return raytracer::Color().R(0.0).G(0.0).B(0.0);
}

raytracer::Color Object::GetProceduralColor(raytracer::Vector position) {
	return raytracer::Color().R(1.0).G(1.0).B(1.0);
}