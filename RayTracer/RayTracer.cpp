#define _USE_MATH_DEFINES

#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>		// must be downloaded 
#include <GL/freeglut.h>	// must be downloaded unless you have an Apple
#endif

#include "Camera.h"
#include "Color.h"
#include "Config.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Diamond.h"
#include "Ellipsoid.h"
#include "FileUtils.h"
#include "HeartShape.h"
#include "Intersection.h"
#include "LightSource.h"
#include "Object.h"
#include "PhotonMap.h"
#include "Program.h"
#include "Polinom.h"
#include "Ray.h"
#include "RealRoots.h"
#include "Rectangle.h"
#include "Scene.h"
#include "Triangle.h"
#include "TriangleMesh.h"
#include "Vector.h"

// IMAGE
float finalImage[KEPSZELESSEG*KEPMAGASSAG * 3];

// raytracer::Colors
raytracer::Color black(0.0, 0.0, 0.0);
raytracer::Color darkGray(0.1f, 0.1f, 0.1f);
raytracer::Color yellow(0.98f, 0.94f, 0.24f);
raytracer::Color encyan(0.31f, 0.91f, 0.91f);
raytracer::Color green(0.0, 0.1f, 0.0);
raytracer::Color red(1.0, 0.0, 0.0);
raytracer::Color white(1.0, 1.0, 1.0);
raytracer::Color blue(0.0, 0.0, 1.0);
raytracer::Color lightBlue(0.41f, 0.79f, 0.91f);
raytracer::Color brown(0.78f, 0.58f, 0.47f);

// n and k values
raytracer::Color copperN(0.2f, 1.1f, 1.2f);
raytracer::Color copperK(3.6f, 2.6f, 2.3f);
raytracer::Color silverN(0.14f, 0.16f, 0.13f);
raytracer::Color silverK(4.1f, 2.3f, 3.1f);
raytracer::Color goldN(0.17f, 0.35f, 1.5f);
raytracer::Color goldK(3.1f, 2.7f, 1.9f);
raytracer::Color transparentN(1.1f, 1.1f, 1.1f);
raytracer::Color transparentK(0.0, 0.0, 0.0);
raytracer::Color glassN(1.5f, 1.5f, 1.5f);
raytracer::Color glassK(0.0, 0.0, 0.0);
raytracer::Color diamondN(2.4f, 2.4f, 2.4f);
raytracer::Color diamondK(0.0, 0.0, 0.0);

// LIGHT1
raytracer::LightSource light1;
//Vector light1Position(800.0, 60.0, 500.0);
raytracer::Vector light1Position(400.0, 150.0, 500.0);
double light1Intensity = 500000.0;
raytracer::Color light1Color = white;

// OBJ1 = TABLE
raytracer::Rectangle table;
raytracer::PhotonMap tablePhotonMap;

raytracer::Vector corner1(-600.0, 300.0);
raytracer::Vector corner2(300.0, -600.0);
raytracer::Vector corner3(1200.0, 300.0);
raytracer::Vector corner4(300.0, 1200.0);

// OBJ2 = Golden ring
raytracer::Cylinder goldenRing;

raytracer::Vector goldenRingCenter(200.0, 100.0, 0.0);
double goldernRingRadius = 300.0;
double goldernRingHeight = 125.0;

// OBJ3 = Silver ring
raytracer::Cylinder silverRing;

raytracer::Vector silverRingCenter(450.0, 100.0, 0.0);
double silverRingRadius = 300.0;
double silverRingHeight = 125.0;

// OBJ4 = Copper ring
raytracer::Cylinder copperRing;

raytracer::Vector copperRingCenter(150.0, 100.0, 0.0);
double copperRingRadius = 300.0;
double copperRingHeight = 100.0;

// OBJ5 = Copper disk
raytracer::Cylinder copperDisk;

raytracer::Vector copperDiskCenter(0.0, 800.0, 0.0);
double copperDiskRadius = 200.0;
double copperDiskHeight = 75.0;

// OBJ6 = Diamond disk
raytracer::Cylinder diamondDisk;

//Vector diamondDiskCenter(550.0, 100.0, 0.0);
raytracer::Vector diamondDiskCenter(200.0, 150.0, 0.0);
double diamondDiskRadius = 150.0;
double diamondDiskHeight = 400.0;

// OBJ7 = Diamond1
raytracer::Diamond diamond1;
raytracer::Ellipsoid diamond1Bounding;
raytracer::TriangleMesh diamond1TriangleMesh;
raytracer::Triangle diamond1Triangles[MAX_HAROMSZOG_SZAM];

raytracer::Vector diamond1Center(550.0, 0.0, 130.0);
double diamond1Width = 125.0;
double diamond1Depth = 125.0;
double diamond1Height = 100.0;

// OBJ8 = Diamond2
raytracer::Diamond diamond2;
raytracer::Ellipsoid diamond2Bounding;
raytracer::TriangleMesh diamond2TriangleMesh;
raytracer::Triangle diamond2Triangles[MAX_HAROMSZOG_SZAM];

raytracer::Vector diamond2Center(0, 100.0, 115.0);
double diamond2Width = 125.0;
double diamond2Depth = 125.0;
double diamond2Height = 100.0;

// OBJ9 = Diamond Ellipsoid
raytracer::Ellipsoid diamondElli;
raytracer::Vector diamondElliCenter(200.0, 150.0, 200.0);
double diamondElliWidth = 75.0;
double diamondElliDepth = 75.0;
double diamondElliHeight = 200.0;

// OBJ10 = Diffuse Ellipsoid
raytracer::Ellipsoid diffuseElli;
raytracer::Vector diffuseElliCenter(200.0, 150.0, 200.0);
double diffuseElliWidth = 75.0;
double diffuseElliDepth = 75.0;
double diffuseElliHeight = 200.0;

// OBJ11 = Heart Shape
raytracer::HeartShape heart;

// CAMERA
raytracer::Camera camera;
//Vector lookAt(300.0,-300.0,300.0);
double r = 1.0;
double angle = 45;
double diffX = r * cos(angle * 180 / PI);
double diffY = r * sin(angle * 180 / PI);
//Vector lookAt(0.0 + diffX, 575.0 + diffY, 0.0);
raytracer::Vector lookAt(0.0, 583.0, 0.0);
double planeWidth = 600.0;
double planeHeight = 600.0;
double planeAngle1 = 0.0;
//double planeAngle2 = 25.0;
double planeAngle2 = 0.0;
double fovDegree = 54.0;
//double fovDegree = 90.0;

// SCENE
raytracer::Scene scene;
raytracer::Color ambient = lightBlue;

void buildDiamondTriangleMesh(raytracer::TriangleMesh *mesh, raytracer::Triangle *triangles, raytracer::Vector center, double width, double depth, double height) {
	// Diamond1

	mesh->Clear();

	int i = 0, u, v;
	double uStart, uEnd, vStart, vEnd;

	raytracer::Vector vertex1, vertex2, vertex3, vertex4;

	double fuggStart = KEZDO_U;
	double fuggInterval = 1.0 - fuggStart;
	double fuggUnit = fuggInterval / TESSZELLACIO_U;

	uStart = KEZDO_U;
	vStart = 0.0;

	vertex1.X(width * sin(uStart*PI) * cos(vStart * 2 * PI));
	vertex1.Y(depth * sin(uStart*PI) * sin(vStart * 2 * PI));
	vertex1.Z(height * cos(uStart*PI));
	vertex1 = center + vertex1;

	for (v = 1; v < TESSZELLACIO_V - 1; v++) {
		vStart = (double)v / TESSZELLACIO_V;
		vEnd = (double)(v + 1) / TESSZELLACIO_V;

		vertex2.X(width * sin(uStart*PI) * cos(vStart * 2 * PI));
		vertex2.Y(depth * sin(uStart*PI) * sin(vStart * 2 * PI));
		vertex2.Z(height * cos(uStart*PI));
		vertex2 = center + vertex2;

		vertex3.X(width * sin(uStart*PI) * cos(vEnd * 2 * PI));
		vertex3.Y(depth * sin(uStart*PI) * sin(vEnd * 2 * PI));
		vertex3.Z(height * cos(uStart*PI));
		vertex3 = center + vertex3;

		triangles[i].Vertex1(vertex1).Vertex2(vertex2).Vertex3(vertex3);
		mesh->AddTriangle(&triangles[i]);

		i++;
	}

	for (u = 0; u < (TESSZELLACIO_U - 1); u++) {
		uStart = fuggStart + pow(u * fuggUnit, TESSZ_HATVANY);
		uEnd = fuggStart + pow((u + 1) * fuggUnit, TESSZ_HATVANY);

		for (v = 0; v < TESSZELLACIO_V; v++) {

			vStart = (double)v / TESSZELLACIO_V;
			vEnd = (double)(v + 1) / TESSZELLACIO_V;

			vertex1.X(width * sin(uStart*PI) * cos(vStart * 2 * PI));
			vertex1.Y(depth * sin(uStart*PI) * sin(vStart * 2 * PI));
			vertex1.Z(height * cos(uStart*PI));
			vertex1 = center + vertex1;

			vertex2.X(width * sin(uStart*PI) * cos(vEnd * 2 * PI));
			vertex2.Y(depth * sin(uStart*PI) * sin(vEnd * 2 * PI));
			vertex2.Z(height * cos(uStart*PI));
			vertex2 = center + vertex2;

			vertex3.X(width * sin(uEnd*PI) * cos(vEnd * 2 * PI));
			vertex3.Y(depth * sin(uEnd*PI) * sin(vEnd * 2 * PI));
			vertex3.Z(height * cos(uEnd*PI));
			vertex3 = center + vertex3;

			vertex4.X(width * sin(uEnd*PI) * cos(vStart * 2 * PI));
			vertex4.Y(depth * sin(uEnd*PI) * sin(vStart * 2 * PI));
			vertex4.Z(height * cos(uEnd*PI));
			vertex4 = center + vertex4;

			triangles[i].Vertex1(vertex1).Vertex2(vertex2).Vertex3(vertex3);
			mesh->AddTriangle(&triangles[i]);

			i++;

			triangles[i].Vertex1(vertex1).Vertex2(vertex3).Vertex3(vertex4);
			mesh->AddTriangle(&triangles[i]);

			i++;
		}
	}

	uStart = fuggStart + pow((TESSZELLACIO_U - 1) * fuggUnit, TESSZ_HATVANY);
	uEnd = 1.0;

	for (v = 0; v < TESSZELLACIO_V; v++) {
		vStart = (double)v / TESSZELLACIO_V;
		vEnd = (double)(v + 1) / TESSZELLACIO_V;

		vertex1.X(width * sin(uStart*PI) * cos(vStart * 2 * PI));
		vertex1.Y(depth * sin(uStart*PI) * sin(vStart * 2 * PI));
		vertex1.Z(height * cos(uStart*PI));
		vertex1 = center + vertex1;

		vertex2.X(width * sin(uStart*PI) * cos(vEnd * 2 * PI));
		vertex2.Y(depth * sin(uStart*PI) * sin(vEnd * 2 * PI));
		vertex2.Z(height * cos(uStart*PI));
		vertex2 = center + vertex2;

		vertex3.X(width * sin(uEnd*PI) * cos(vStart * 2 * PI));
		vertex3.Y(depth * sin(uEnd*PI) * sin(vStart * 2 * PI));
		vertex3.Z(height * cos(uEnd*PI));
		vertex3 = center + vertex3;

		triangles[i].Vertex1(vertex1).Vertex2(vertex2).Vertex3(vertex3);
		mesh->AddTriangle(&triangles[i]);

		i++;
	}

}

void initLights() {
	// LIGHT1
	light1.LightPosition(light1Position).LightColor(light1Color).Intensity(light1Intensity);
}

void initObjects() {
	// OBJ1 = ASZTAL
	//asztal.Diffuse().Kd(darkGray).Ka(darkGray * PI).Ks(white).Shine(1000.0).ProcMode(raytracer::Object::ProceduralMode::MULTIPLY);
	table.Diffuse().Kd(darkGray).Ka(darkGray * PI).Ks(white).Shine(1000.0).ProcMode(raytracer::Object::ProceduralMode::NONE);
	table.Corner1(corner1).Corner2(corner2).Corner3(corner3).Corner4(corner4);
	table.SetPhotonMap(&tablePhotonMap);

	// OBJ2 = Golden ring
	goldenRing.N(goldN).K(goldK).Reflective();
	goldenRing.BasePoint(goldenRingCenter).Radius(goldernRingRadius).Height(goldernRingHeight);

	// OBJ3 = Silver ring
	silverRing.N(silverN).K(silverK).Reflective();
	silverRing.BasePoint(silverRingCenter).Radius(silverRingRadius).Height(silverRingHeight);

	// OBJ4 = Copper ring
	copperRing.N(copperN).K(copperK).Reflective();
	copperRing.BasePoint(copperRingCenter).Radius(copperRingRadius).Height(copperRingHeight);

	// OBJ5 = Copper disk
	copperDisk.N(copperN).K(copperK).Reflective();
	copperDisk.BasePoint(copperDiskCenter).Radius(copperDiskRadius).Height(copperDiskHeight).Solid();

	// OBJ6 = Diamond disk
	diamondDisk.N(diamondN).K(diamondK).Reflective().Refractive();
	diamondDisk.BasePoint(diamondDiskCenter).Radius(diamondDiskRadius).Height(diamondDiskHeight).Solid();

	// OBJ7 = Diamond1
	buildDiamondTriangleMesh(&diamond1TriangleMesh, diamond1Triangles, diamond1Center, diamond1Width, diamond1Depth, diamond1Height);
	diamond1Bounding.Center(diamond1Center).A(diamond1Width).B(diamond1Depth).C(diamond1Height);
	diamond1.N(diamondN).K(diamondK).Reflective().Refractive();
	diamond1.BoundingVolume(&diamond1Bounding).DiamondBody(&diamond1TriangleMesh);

	// OBJ8 = Diamond2
	buildDiamondTriangleMesh(&diamond2TriangleMesh, diamond2Triangles, diamond2Center, diamond2Width, diamond2Depth, diamond2Height);
	diamond2Bounding.Center(diamond2Center).A(diamond2Width).B(diamond2Depth).C(diamond2Height);
	diamond2.N(diamondN).K(diamondK).Reflective().Refractive();
	diamond2.BoundingVolume(&diamond2Bounding).DiamondBody(&diamond2TriangleMesh);

	// OBJ9 = Diamond Ellipsoid
	//diamondElli.N(diamondN).K(diamondK).Reflective().Refractive();
	diamondElli.N(transparentN).K(transparentK).Reflective().Refractive();
	diamondElli.Center(diamondElliCenter).A(diamondElliWidth).B(diamondElliDepth).C(diamondElliHeight);

	// OBJ10 = Diffuse Ellipsoid
	//diffuseElli.Kd(green).Ka(green * PI).Diffuse();
	diffuseElli.Kd(green).Ka(black).Diffuse();
	diffuseElli.Center(diffuseElliCenter).A(diffuseElliWidth).B(diffuseElliDepth).C(diffuseElliHeight);

	heart.Kd(lightBlue).Ka(lightBlue / PI).Diffuse();
	//heart.N(diamondN).K(diamondK).Reflective().Refractive();

}

void initCamera() {
	double cos1 = cos((planeAngle1 / 180.0) * PI);
	double sin1 = sin((planeAngle1 / 180.0) * PI);
	double cos2 = cos((planeAngle2 / 180.0) * PI);
	double sin2 = sin((planeAngle2 / 180.0) * PI);

	double halfPlaneWidth = planeWidth / 2;
	double halfPlaneHeight = planeHeight / 2;
	raytracer::Vector right(halfPlaneWidth * cos1, halfPlaneWidth * sin1, 0.0);
	raytracer::Vector up(0.0, halfPlaneHeight * sin2, halfPlaneHeight * cos2);

	camera.LookAt(lookAt).Right(right).Up(up).FovDegree(fovDegree);
}

void initScene() {
	scene.Cam(&camera).Ambient(ambient);

	scene.AddLightSource(&light1);

	//scene.AddObject(&table);
	//scene.AddObject(&goldenRing);
	//scene.AddObject(&silverRing);
	//scene.AddObject(&copperRing);
	//scene.AddObject(&copperDisk);

	//scene.AddObject(&diamondDisk);

	//scene.AddObject(&diamond1);
	//scene.AddObject(&diamond2);
	//scene.AddObject(&diamondElli);
	//scene.AddObject(&diffuseElli);

	scene.AddObject(&heart);
}

int main(int argc, char **argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(KEPSZELESSEG, KEPMAGASSAG);
	glutInitWindowPosition(100, 100);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

	glutCreateWindow("Ray tracer");

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	initialize();

	glutDisplayFunc(onDisplay);
	glutMouseFunc(onMouseEvent);
	glutKeyboardFunc(onKeyboardEvent);

	glutMainLoop();

	return 0;
}

void initialize() {
	initLights();
	initObjects();
	initCamera();
	initScene();

	scene.Render(finalImage);
	if (TONEMAP) {
		scene.ToneMap(finalImage);
	}

	if (FILE_IRAS) {
		raytracer::files::saveToBmp(FILE_NEV, finalImage, KEPSZELESSEG, KEPMAGASSAG);
	}

	glutPostRedisplay();
}

void onDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(KEPSZELESSEG, KEPMAGASSAG, GL_RGB, GL_FLOAT, finalImage);
	glutSwapBuffers();
}

void onKeyboardEvent(unsigned char key, int x, int y) {
	if (key == 'd') {
		/*
		scene.Render(finalImage);
		scene.ToneMap(finalImage);
		*/
		glutPostRedisplay();
	}
	else if (key == ' ') {
		exit(0);
	}
}

void onMouseEvent(int button, int state, int x, int y) {
	//if (button == GLUT_LEFT && state == GLUT_DOWN);  // A GLUT_LEFT_BUTTON / GLUT_RIGHT_BUTTON illetve GLUT_DOWN / GLUT_UP
}