#ifndef LIGHT_SOURCE_H
#define LIGHT_SOURCE_H

#include "Color.h"
#include "Vector.h"

namespace raytracer {

	class LightSource {
		raytracer::Color lightColor;
		double intensity;
		raytracer::Vector lightPosition;
	public:
		LightSource(double newIntensity = 0.0, raytracer::Vector newPosition = raytracer::Vector().X(0.0).Y(0.0).Z(0.0), raytracer::Color newColor = raytracer::Color().R(1.0).G(1.0).B(1.0));

		LightSource& LightColor(raytracer::Color newColor);

		raytracer::Color GetLightColor();

		LightSource& Intensity(double newIntensity);

		double GetIntensity();

		LightSource& LightPosition(raytracer::Vector newPosition);

		raytracer::Vector getLightPosition();

		raytracer::Color GetLightAt(raytracer::Vector point);

	};

}

#endif // !LIGHT_SOURCE_H