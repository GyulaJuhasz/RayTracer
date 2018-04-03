#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>
#include <GL/freeglut.h>
#endif

#include <iostream>

#include "Scene.h"

using raytracer::Scene;

Scene::Scene() {
	ambient.R(0.0).G(0.0).B(0.0);
	camera = NULL;
	ClearObjects();
	ClearLights();
}

Scene& Scene::Ambient(raytracer::Color newAmbient) {
	ambient = newAmbient;
	return *this;
}

raytracer::Color Scene::GetAmbient() { return ambient; }

Scene& Scene::Cam(raytracer::Camera *newCamera) {
	camera = newCamera;
	return *this;
}

raytracer::Camera* Scene::GetCam() { return camera; }

Scene& Scene::AddObject(raytracer::Object *newObject) {
	if (numOfObjects < MAX_OBJ_SZAM) {
		objects[numOfObjects++] = newObject;
	}
	return *this;
}

int Scene::GetNumberOfObjects() { return numOfObjects; }

void Scene::ClearObjects() {
	for (int i = 0; i < MAX_OBJ_SZAM; i++) objects[i] = NULL;
	numOfObjects = 0;
}

Scene& Scene::AddLightSource(raytracer::LightSource *newLight) {
	if (numOfLights < MAX_FENY_SZAM) {
		lights[numOfLights++] = newLight;
	}
	return *this;
}

int Scene::GetNumberOfLightSources() { return numOfLights; }

void Scene::ClearLights() {
	for (int i = 0; i < MAX_FENY_SZAM; i++) lights[i] = NULL;
	numOfLights = 0;
}

raytracer::Intersection Scene::GetFirstIntersect(raytracer::Ray ray) {
	raytracer::Intersection result;
	raytracer::Intersection actual;
	double tmin = 99999999999.0;
	if (numOfObjects > 0) {
		raytracer::Intersection temp;
		for (int i = 0; i < numOfObjects; i++) {
			actual = objects[i]->Intersect(ray);
			if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
				tmin = actual.GetParam();
				result = actual;
				result.Index(i);
			}
		}
	}

	return result;
}

raytracer::Color Scene::GetDirectLight(raytracer::Vector point, raytracer::Vector N, raytracer::Vector V, int index) {
	// point = point on surface
	// N = normal on surface
	// V = view direction = vector from point to the eye of the camera
	// index = index of the object thats lighened raytracer::Color we are calculating

	// outerPoint = point + epsilon * normal
	raytracer::Vector outerPoint = point + N * EPSZILON_TOLAS;

	raytracer::Object *object = objects[index];

	raytracer::Color ka = object->GetKa();
	raytracer::Color kd = object->GetKd();
	raytracer::Color ks = object->GetKs();
	double shine = object->GetShine();
	raytracer::Color procedural = object->GetProceduralColor(point);
	raytracer::Object::ProceduralMode proceduralMode = object->GetProcMode();

	switch (proceduralMode)
	{
	case raytracer::Object::ProceduralMode::ADD:
		kd = kd + procedural;
		ka = ka + procedural * PI;
		break;
	case raytracer::Object::ProceduralMode::MULTIPLY:
		kd = kd * procedural;
		ka = ka * procedural * PI;
		break;
	case raytracer::Object::ProceduralMode::REPLACE:
		kd = procedural;
		ka = procedural * PI;
		break;
	case raytracer::Object::ProceduralMode::NONE:
	default:
		break;
	}

	// result = Ka * La
	raytracer::Color result = ka * ambient;

	raytracer::Vector lightPosition;
	// L = direction of shadow ray = vector from outerPoint to actual lightsource
	raytracer::Vector L;
	double lightDistance;
	raytracer::Intersection lightIntersection;
	raytracer::Vector intersectVector;

	raytracer::LightSource *actualLightSource;

	for (int i = 0; i < numOfLights; i++) {
		actualLightSource = lights[i];
		lightPosition = actualLightSource->getLightPosition();
		L = lightPosition - outerPoint;
		lightDistance = L.GetLenght();

		raytracer::Ray shadowRay(outerPoint, L);
		lightIntersection = GetFirstIntersect(shadowRay);
		intersectVector = lightIntersection.GetIntersectPoint() - point;

		if ((lightIntersection.GetIndex() == -1) || (intersectVector.GetLenght() > lightDistance)) {
			L.Normalize();
			raytracer::Vector H = (L + V);
			H.Normalize();
			double costheta = N * L;
			double cosdelta = N * H;

			// Lin: Light Source's raytracer::Color in point
			raytracer::Color Lin = actualLightSource->GetLightAt(point);

			if (costheta > 0.0) {
				// Result = Ka * La + kd * Lin * costetha
				result = result + kd * (Lin * costheta);

				if (cosdelta > 0.0) {
					// Result = Ka * La + kd * Lin * costetha + ks * Lin * cosdelta^shine
					result = result + ks * (Lin * pow(cosdelta, shine));
				}
			}
		}
	}

	// Adding caustics
	if (KAUSZTIKA) result = result + (object->GetPhoton(point) * KAUSZTIKA_SZORZO);

	return result;
}

raytracer::Color Scene::Trace(raytracer::Ray ray, int depth, bool inside) {
	if (depth > D_MAX) return ambient;

	raytracer::Intersection firstIntersection = GetFirstIntersect(ray);
	int index = firstIntersection.GetIndex();
	if (index == -1) return ambient;

	raytracer::Object *object = objects[index];

	raytracer::Vector N = firstIntersection.GetNormal();
	// IntersectPoint = x
	raytracer::Vector point = firstIntersection.GetIntersectPoint();
	raytracer::Vector V = (camera->GetEye() - point).GetNormalized();
	// RAy DIR = v
	raytracer::Vector rayDir = ray.GetV().GetNormalized();
	if (N * rayDir > 0) N = N * -1.0;

	// Result = 0
	raytracer::Color result;

	// If diffuse: Result = Result + DirectLight
	if (object->IsDiffuse() && LAMPAFENY) {
		bool LOL;
		if (depth > 0) {
			LOL = true;
		}
		result = result + GetDirectLight(point, N, V, index);
	}

	raytracer::Color objectN;
	raytracer::Color fresnel;
	raytracer::Color kt;
	double n;
	if (object->IsReflective() || object->IsRefractive()) {
		fresnel = object->GetFresnel(N, rayDir);
		kt.R(1.0 - fresnel.GetR()).G(1.0 - fresnel.GetG()).B(1.0 - fresnel.GetB());
		objectN = object->GetN();
		n = objectN.GetR();
		if (inside) {
			objectN.R(1.0 / objectN.GetR()).G(1.0 / objectN.GetG()).B(1.0 / objectN.GetB());
			n = 1.0 / n;
		}
	}
	if (object->IsRefractive() && TORES) {
		bool fullReflect = false;
		if (inside) {
			double cosalfa = N * (rayDir * (-1.0));
			double sinalfanegyzet = 1.0 - cosalfa*cosalfa;
			double cosbetanegyzet = 1.0 - sinalfanegyzet / (n*n);
			if (cosbetanegyzet < 0.0) {
				fullReflect = true;
				fresnel.R(1.0).G(1.0).B(1.0);
			}
		}

		if (!fullReflect) {
			// If refractive: Result = Result + Refracted * kt
			result = result + Trace(ray.Refract(point, N, n), depth + 1, !inside) * kt;
		}
	}

	if (object->IsReflective() && TUKROZES) {
		// If refractive: Result = Result + Reflected * Fresnel
		result = result + Trace(ray.Reflect(point, N), depth + 1, inside) * fresnel;
	}

	return result;
}

void Scene::Shoot(raytracer::Color intensity, raytracer::Ray ray, int depth, bool inside) {
	if (depth > D_MAX) return;

	raytracer::Intersection firstIntersection = GetFirstIntersect(ray);
	int index = firstIntersection.GetIndex();
	if (index == -1) return;

	raytracer::Object *object = objects[index];

	raytracer::Vector N = firstIntersection.GetNormal();
	// IntersectPoint = x
	raytracer::Vector point = firstIntersection.GetIntersectPoint();
	// RAy DIR = v
	raytracer::Vector rayDir = ray.GetV().GetNormalized();
	if (N * rayDir > 0) N = N * -1.0;

	raytracer::Color objectN;
	raytracer::Color fresnel;
	raytracer::Color kt;
	double n;
	if (object->IsReflective() || object->IsRefractive()) {
		fresnel = object->GetFresnel(N, rayDir);
		kt.R(1.0 - fresnel.GetR()).G(1.0 - fresnel.GetG()).B(1.0 - fresnel.GetB());
		objectN = object->GetN();
		n = objectN.GetR();
		if (inside) {
			objectN.R(1.0 / objectN.GetR()).G(1.0 / objectN.GetG()).B(1.0 / objectN.GetB());
			n = 1.0 / n;
		}
	}
	if (object->IsRefractive() && TORES) {
		bool fullReflect = false;
		if (inside) {
			double cosalfa = N * (rayDir * (-1.0));
			double sinalfanegyzet = 1.0 - cosalfa*cosalfa;
			double cosbetanegyzet = 1.0 - sinalfanegyzet / (n*n);
			if (cosbetanegyzet < 0.0) {
				fullReflect = true;
				fresnel.R(1.0).G(1.0).B(1.0);
			}
		}

		if (!fullReflect) {
			// If refractive: Shoot (intensity * kt) in refract dir
			Shoot(intensity * kt, ray.Refract(point, N, n), depth + 1, !inside);
		}
	}

	if (object->IsReflective() && TUKROZES) {
		// If reflective: Shoot (intensity * fresnel) in reflect dir
		Shoot(intensity * fresnel, ray.Reflect(point, N), depth + 1, inside);
	}

	// If diffuse, and photon has hit at least one other object, store it in photonMap
	if (object->IsDiffuse() && LAMPAFENY && depth > 0) object->AddPhoton(intensity, point);
	return;
}

void Scene::ToneMap(float *image) {

	double sumOfLuminance = 0.0;
	bool needed = false;

	for (int i = 0; i < KEPSZELESSEG*KEPMAGASSAG * 3; i++) {
		if (image[i] > 1.0) {
			needed = true;
			break;
		}
	}

	if (!needed) return;

	for (int i = 0; i < KEPSZELESSEG*KEPMAGASSAG * 3; i += 3) {
		sumOfLuminance += 0.21f * image[i];
		sumOfLuminance += 0.72f * image[i + 1];
		sumOfLuminance += 0.07f * image[i + 2];
	}

	double averageOfLuminance = sumOfLuminance / (KEPSZELESSEG*KEPMAGASSAG);
	double I, Ir, D;

	for (int i = 0; i < KEPSZELESSEG*KEPMAGASSAG * 3; i += 3) {
		I = 0.21f*image[i] + 0.72f*image[i + 1] + 0.07f*image[i + 2];
		if (I < 0.0) {
			image[i] = 0.0;
			image[i + 1] = 0.0;
			image[i + 2] = 0.0;
		}
		else {
			Ir = I / averageOfLuminance;
			D = (ALFA*Ir) / (1 + ALFA*Ir);
			image[i] = (float)(image[i] * D / I);
			image[i + 1] = (float)(image[i + 1] * D / I);
			image[i + 2] = (float)(image[i + 2] * D / I);
		}
	}
}

void Scene::Render(float *image) {
	raytracer::Ray ray;
	if (camera != NULL) {
		if (KAUSZTIKA) {
			if (FOTON_FAJLBOL) {

				std::cout << std::endl << std::endl << "Fotonterkepek fajlbol..." << std::endl << std::endl;

				for (int i = 0; i < numOfObjects; i++) {
					if (objects[i]->IsDiffuse()) {
						char fileName[200];
						sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", FOTONSZAM, i);
						//sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", fotonszam, i);
						objects[i]->ReadPhotonMapFromFile(fileName);
					}
				}

				std::cout << "Kesz." << std::endl << std::endl;

			}
			else {

				std::cout << "Fotonloves start..." << std::endl << std::endl;

				double fx, fy, fz;
				double length, intensity;
				raytracer::Ray photonRay;
				raytracer::Color lightColor;
				raytracer::Vector lightPosition;
				raytracer::Vector direction;

				raytracer::Intersection testIntersection;

				raytracer::LightSource *actual;

				for (int i = 0; i < numOfLights; i++) {
					std::cout << i + 1 << ". lampa fotonloves..." << std::endl;

					actual = lights[i];

					lightPosition = actual->getLightPosition();
					intensity = actual->GetIntensity() / FOTONSZAM;
					lightColor = actual->GetLightColor();

					//intenzitas/=fotonszam;
					//intenzitas*=KAUSZTIKA_SZORZO;

					bool good = false;

					unsigned long starttime = glutGet(GLUT_ELAPSED_TIME);		// program inditasa ota eltelt ido
					unsigned long actualtime;
					unsigned long ellapsed;
					unsigned long remaining;
					unsigned long percent;

					for (int j = 0; j < FOTONSZAM; j++) {
						//for(unsigned long j=0;j<fotonszam;j++) {
						//if ( !(j% (FOTONSZAM / 100 )) ) cout << (double)(j*100) / FOTONSZAM << "%..." << endl;

						if (!(j % (FOTONSZAM / 100))) {
							//cout << j << endl;
							percent = (unsigned long)(j * 100) / FOTONSZAM;
							std::cout << percent << "%";
							if (percent > 0) {
								actualtime = glutGet(GLUT_ELAPSED_TIME);		// program inditasa ota eltelt ido
								ellapsed = actualtime - starttime;
								remaining = (ellapsed / percent) * (100 - percent);
								std::cout << ", becsult hatralevo ido: " << remaining / 1000 << " s" << std::endl;
							}
							else std::cout << std::endl;
						}
						do {

							fx = RND;
							fy = RND;
							fz = RND;
							if (fz > 0.0) continue;

							if ((fx*fx + fy*fy + fz*fz) > 1.0) continue;

							length = sqrt(fx*fx + fy*fy + fz*fz);

							if (length > 0.0) {
								fx /= length;
								fy /= length;
								fz /= length;
							}

							direction.X(fx).Y(fy).Z(fz);
							photonRay.P0(lightPosition).V(direction);
							testIntersection = GetFirstIntersect(photonRay);
							if (testIntersection.GetIndex() == -1) continue;
							good = true;
						} while (!good);
						Shoot(lightColor * intensity, photonRay, 0, false);
					}
				}

				std::cout << std::endl << std::endl << "Fotonloves kesz." << std::endl << std::endl;

				if (FOTON_FAJLBA) {

					std::cout << std::endl << std::endl << "Fotonterkepek fajlba..." << std::endl << std::endl;

					for (int i = 0; i < numOfObjects; i++) {
						if (objects[i]->IsDiffuse()) {
							char fileName[200];
							sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", FOTONSZAM, i);
							//sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", fotonszam, i);
							objects[i]->WritePhotonMapToFile(fileName);
						}
					}
					std::cout << "Kesz." << std::endl << std::endl;
				}

			}
		}
		std::cout << "Sugarkovetes start..." << std::endl << std::endl;

		raytracer::Color actual;
		for (int z = 0; z < KEPMAGASSAG; z++) {

			std::cout << (z + 1) << ". sor..." << std::endl;

			for (int x = 0; x < KEPSZELESSEG; x++) {
				bool willFail = false;
				ray = camera->GetRay(x, z);

				if (x == 260 && z == 175) {
					willFail = true;
				}

				actual = Trace(ray, 0, false);

				int actualIndex = (KEPMAGASSAG - 1 - z) * 3 * KEPSZELESSEG + x * 3;
				image[actualIndex] = (float)actual.GetR();
				image[actualIndex + 1] = (float)actual.GetG();
				image[actualIndex + 2] = (float)actual.GetB();
			}
		}
	}
}