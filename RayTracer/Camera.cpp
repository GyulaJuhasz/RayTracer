#include <cmath>

#include "Camera.h"
#include "Config.h"

using raytracer::Camera;

Camera::Camera() {
	lookAtPoint.X(0.0).Y(0.0).Z(0.0).H(1.0);
	right.X(0.0).Y(0.0).Z(0.0).H(1.0);
	up.X(0.0).Y(0.0).Z(0.0).H(1.0);
	eye.X(0.0).Y(0.0).Z(0.0).H(1.0);
	fovDegree = 0.0;
}

Camera& Camera::LookAt(raytracer::Vector newLookAt) {
	lookAtPoint = newLookAt;
	CalculateEye();
	return *this;
}

raytracer::Vector Camera::GetLookAtPoint() { return lookAtPoint; }

Camera& Camera::Right(raytracer::Vector newRight) {
	right = newRight;
	CalculateEye();
	return *this;
}

raytracer::Vector Camera::GetRightVector() { return right; }

Camera& Camera::Up(raytracer::Vector newUp) {
	up = newUp;
	CalculateEye();
	return *this;
}

raytracer::Vector Camera::GetupVector() { return up; }

Camera& Camera::FovDegree(double newFov) {
	fovDegree = newFov;
	CalculateEye();
	return *this;
}

double Camera::GetFovDegree() { return fovDegree; }

raytracer::Vector Camera::GetEye() { return eye; }

raytracer::Ray Camera::GetRay(int x, int z) {
	raytracer::Ray result;
	result.P0(eye).V(CenterOfPixel(x, z) - eye);
	return result;
}

void Camera::CalculateEye() {
	double halfFovRad = ((fovDegree / 2.0) / 180.0) * PI;
	double tanHalfFov = tan(halfFovRad);
	double eyeDistance = (double)KEPMAGASSAG / (2.0 * tanHalfFov);
	raytracer::Vector cameraAxisNormal = (right % up).GetNormalized();
	eye = lookAtPoint + cameraAxisNormal * eyeDistance;
}

raytracer::Vector Camera::CenterOfPixel(int x, int z) {
	raytracer::Vector result;
	//result = lookAtPoint + right * ( ((double)2*x/KEPSZELESSEG) - 1.0 ) + up * ( ((double)2*z/KEPMAGASSAG)-1 );
	result = lookAtPoint + right * (((double)2 * x / KEPSZELESSEG) - 1.0) + up * (1.0 - ((double)(2 * z - 2) / KEPMAGASSAG));
	return result;
}