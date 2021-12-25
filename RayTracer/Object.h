#ifndef OBJECT_H
#define OBJECT_H

#include "Color.h"
#include "Intersection.h"
#include "PhotonMap.h"
#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Object {

	public:
		enum ProceduralMode {
			NONE,
			ADD,
			MULTIPLY,
			REPLACE
		};

	protected:
		bool isDiffuse, isReflective, isRefractive;
		raytracer::Color kd, ks, ka, F0, n, k;
		double shine;
		ProceduralMode proceduralMode;
		raytracer::PhotonMap* photonMap;
		Object *boundingObject;
		virtual raytracer::Intersection SpecificIntersect(raytracer::Ray ray) = 0;

	public:

		Object();

		Object& Diffuse(bool diffuse = true);

		bool IsDiffuse();

		Object& Reflective(bool reflective = true);

		bool IsReflective();

		Object& Refractive(bool refractive = true);

		bool IsRefractive();

		Object& Kd(raytracer::Color newKd);

		raytracer::Color GetKd();

		Object& Ks(raytracer::Color newKs);

		raytracer::Color GetKs();

		Object& Ka(raytracer::Color newKa);

		raytracer::Color GetKa();

		Object& Shine(double newShine);

		double GetShine();

		void calculateFresnel();

		Object& N(raytracer::Color newN);

		raytracer::Color GetN();

		Object& K(raytracer::Color newK);

		Object& BoundingObject(raytracer::Object *newBoundingObject);

		raytracer::Color GetK();

		raytracer::Color GetFresnel(raytracer::Vector N, raytracer::Vector V);

		Object& ProcMode(ProceduralMode newMode);

		ProceduralMode GetProcMode();

		raytracer::Intersection Intersect(raytracer::Ray ray);

		void SetPhotonMap(PhotonMap *newMap);

		raytracer::PhotonMap* GetPhotonMap();

		void WritePhotonMapToFile(const char* fileName);

		void ReadPhotonMapFromFile(const char* fileName);

		virtual void AddPhoton(raytracer::Color intensity, raytracer::Vector position);

		virtual raytracer::Color GetPhoton(raytracer::Vector position);

		virtual raytracer::Color GetProceduralColor(raytracer::Vector position);
	};

}

#endif // !OBJECT_H