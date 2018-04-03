#include "LightSource.h"
#include "MathUtils.h"

using raytracer::LightSource;

LightSource::LightSource(double newIntensity, raytracer::Vector newPosition, raytracer::Color newColor):intensity(newIntensity), lightPosition(newPosition), lightColor(newColor) {}

LightSource& LightSource::LightColor(raytracer::Color newColor) {
	lightColor = newColor;
	return *this;
}

raytracer::Color LightSource::GetLightColor() { return lightColor; }

LightSource& LightSource::Intensity(double newIntensity) {
	intensity = newIntensity;
	return *this;
}

double LightSource::GetIntensity() { return intensity; }

LightSource& LightSource::LightPosition(raytracer::Vector newPosition) {
	lightPosition = newPosition;
	return *this;
}

raytracer::Vector LightSource::getLightPosition() { return lightPosition; }

raytracer::Color LightSource::GetLightAt(raytracer::Vector point) {
	raytracer::Vector distanceVector = lightPosition - point;
	double distance = distanceVector.GetLenght();
	double power;
	if (raytracer::math::IsZero(distance)) power = intensity;
	else power = intensity / (distance * distance);
	return lightColor * power;
}
