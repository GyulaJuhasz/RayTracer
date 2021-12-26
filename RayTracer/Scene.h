#ifndef SCENE_H
#define SCENE_H

#include "Camera.h"
#include "Color.h"
#include "Config.h"
#include "Intersection.h"
#include "LightSource.h"
#include "Object.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Scene {
		raytracer::Color ambient;
		raytracer::Camera* camera;
		raytracer::Object* objects[MAX_OBJ_SZAM];
		int numOfObjects;
		raytracer::LightSource* lights[MAX_FENY_SZAM];
		int numOfLights;

	public:
		Scene();

		Scene& Ambient(raytracer::Color newAmbient);

		raytracer::Color GetAmbient();

		Scene& Cam(raytracer::Camera *newCamera);

		raytracer::Camera* GetCam();

		Scene& AddObject(raytracer::Object *newObject);

		int GetNumberOfObjects();

		void ClearObjects();

		Scene& AddLightSource(raytracer::LightSource *newLight);

		int GetNumberOfLightSources();

		void ClearLights();

		raytracer::Intersection GetFirstIntersect(raytracer::Ray ray);

		raytracer::Color GetDirectLight(raytracer::Vector point, raytracer::Vector N, raytracer::Vector V, int index);

		raytracer::Color Trace(raytracer::Ray ray, int depth, bool inside);

		void Shoot(raytracer::Color intensity, raytracer::Ray ray, int depth, bool inside);

		void ToneMap(float *image, int imageWidth, int imageHeight);

		void Render(float *image, int imageWidth, int imageHeight);
	};

}

#endif // !SCENE_H