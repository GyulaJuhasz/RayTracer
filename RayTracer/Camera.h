#ifndef CAMERA_H
#define CAMERA_H

#include "Ray.h"
#include "Vector.h"

namespace raytracer {

	class Camera {
		raytracer::Vector lookAtPoint;
		raytracer::Vector right;
		raytracer::Vector up;
		raytracer::Vector eye;
		double fovDegree;

	public:
		Camera();

		Camera& LookAt(raytracer::Vector newLookAt);

		raytracer::Vector GetLookAtPoint();

		Camera& Right(raytracer::Vector newRight);

		raytracer::Vector GetRightVector();

		Camera& Up(raytracer::Vector newUp);

		raytracer::Vector GetupVector();

		Camera& FovDegree(double newFov);

		double GetFovDegree();

		raytracer::Vector GetEye();

		raytracer::Ray GetRay(int x, int z);

	private:
		void CalculateEye();

		raytracer::Vector CenterOfPixel(int x, int z);

	};

}


#endif // !CAMERA_H