#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

#include <iostream>
#include <fstream>

#define new new_nelkul_is_meg_lehet_csinalni

#define EPSZILON 0.000001f
#define EPSZILON_TOLAS 0.1f
#define PI 3.141592f
#define KEPSZELESSEG 600
#define KEPMAGASSAG 600
#define MAX_OBJ_SZAM 100
#define MAX_FENY_SZAM 10
#define ALFA 2.0
#define D_MAX 50
#define TEXTURA true
#define TEXTURA_X 300
#define TEXTURA_Y 300
//#define FOTONTERKEP_X 1000
//#define FOTONTERKEP_Y 1000
#define FOTONTERKEP_X 1000
#define FOTONTERKEP_Y 1000
#define FOTONSZAM 1000000
#define TUKROZES true
#define TORES true
#define LAMPAFENY true
#define KAUSZTIKA false
#define TONEMAP true
#define KAUSZTIKA_SZORZO 0.1f
#define RND (2.0*(double)rand()/RAND_MAX-1.0)
#define TESSZELLACIO_U 5
#define TESSZELLACIO_V 6
#define KEZDO_U 0.2f
#define TESSZ_HATVANY 4.0
#define MAX_HAROMSZOG_SZAM 800
#define MAX_ATLAG 5
#define ATLAG_SUGAR 3
#define FILE_IRAS true
#define FILE_NEV "test.bmp"
#define FOTON_FAJLBOL false
#define FOTON_FAJLBA false

#define MAX_POLINOM_DEGREE 8

using namespace std;

class FloatHelper {
	static const double epsilon;
public:
	static bool Equals(double a, double b) {
		return fabs(a - b) < epsilon;
	}

	static bool NotEquals(double a, double b) {
		return !Equals(a, b);
	}

	static bool IsZero(double a) {
		return Equals(a, 0.0);
	}

	static int Sig(double a) {
		if (a > 0.0000001) {
			return 1;
		}
		else if (a < -0.0000001) {
			return -1;
		}
		else {
			return 0;
		}
	}

};

const double FloatHelper::epsilon = 0.000001;

class RealRoots {
	double realRoots[MAX_POLINOM_DEGREE];
	int numOfRoots;

public:
	RealRoots() {
		numOfRoots = 0;
	}

	void AddRoot(double x0) {
		if (numOfRoots < MAX_POLINOM_DEGREE) {
			realRoots[numOfRoots++] = x0;
		}
	}

	int GetNumberOfRoots() { return numOfRoots; }

	double GetRoot(int index) {
		if (index < 0 || index > numOfRoots) {
			throw "Wrong INDEX!";
		}

		return realRoots[index];
	}

	double GetMinPositiveRootIfAny() {
		if (numOfRoots == 0) return -1.0;

		bool found = false;
		double minPos = 99999999.0;
		double actualRoot;

		for (int i = 0; i < numOfRoots; i++) {
			actualRoot = realRoots[i];
			if (actualRoot > 0.0) {
				found = true;
			}

			if (actualRoot < minPos) {
				minPos = actualRoot;
			}
		}

		return (found ? minPos : -1.0);

	}

	RealRoots& operator=(const RealRoots& r) {
		if (this != &r) {
			numOfRoots = r.numOfRoots;
			for (int i = 0; i <= numOfRoots; i++) {
				realRoots[i] = r.realRoots[i];
			}
		}

		return *this;
	}

	friend ostream& operator<<(ostream& os, const RealRoots& r);


};

ostream& operator<<(ostream& os, const RealRoots& r) {
	if (r.numOfRoots == 0) {
		os << "NO REAL ROOTS!" << std::endl;
	}
	else {
		for (int i = 0; i < r.numOfRoots; i++) {
			os << "x" << (i + 1) << " = " << (r.realRoots[i]) << std::endl;
		}
	}

	return os;
}


class Polinom {
	double coefficients[MAX_POLINOM_DEGREE];
	int degree;

	static const double numeric_epsilon;
	static const int max_wrong_tries;
	static const int min_good_tries;

public:
	Polinom(int newDegree) {
		if (newDegree < 0 || newDegree > MAX_POLINOM_DEGREE) {
			throw "Invalid degree!";
		}

		degree = newDegree;
	}

	Polinom(const Polinom& p) {
		degree = p.degree;
		for (int i = 0; i <= degree; i++) {
			coefficients[i] = p.coefficients[i];
		}
	}

	int GetDegree() { return degree; }

	double GetCoefficient(int whichDegree) const {
		if (whichDegree < 0 || whichDegree > degree) {
			throw "Invalid degree!";
		}

		return coefficients[whichDegree];
	}

	void SetCoeff(int whichDegree, double coeff) {
		if (whichDegree < 0 || whichDegree > degree) {
			throw "Invalid degree!";
		}

		coefficients[whichDegree] = coeff;
	}

	Polinom& SetCoefficient(int whichDegree, double coeff) {
		if (whichDegree < 0 || whichDegree > degree) {
			throw "Invalid degree!";
		}

		coefficients[whichDegree] = coeff;
		return *this;
	}

	Polinom divideByRoot(double x0) {
		if (degree == 0) {
			throw "Cannot divide 0 degree!";
		}

		Polinom result(degree - 1);

		for (int i = degree - 1; i >= 0; i--) {
			double ai = 0.0;

			for (int j = 0; j <= degree - 1 - i; j++) {
				double power;
				power = 1.0;

				for (int k = 0; k < (degree - 1 - i - j); k++) power *= x0;

				ai += GetCoefficient(degree - j) * power;
			}

			result.SetCoefficient(i, ai);
		}

		return result;
	}

	Polinom GetDerived() {
		if (degree == 0) {
			Polinom result(0);
			result.SetCoefficient(0, 0.0);
			return result;
		}

		Polinom result(degree - 1);

		for (int i = degree - 1; i >= 0; i--) {
			result.SetCoefficient(i, GetCoefficient(i + 1) * (i + 1));
		}

		return result;
	}

	Polinom& operator=(const Polinom& p) {
		if (this != &p) {
			degree = p.degree;
			for (int i = 0; i <= degree; i++) {
				coefficients[i] = p.coefficients[i];
			}
		}

		return *this;
	}

	double GetValueAt(double x0) {
		double result = 0.0;
		double power = 1.0;

		for (int i = 0; i <= degree; i++) {
			if (i != 0) power *= x0;
			result += coefficients[i] * power;
		}

		return result;
	}

	RealRoots getRootsWithNewton(double guess) {
		RealRoots result;

		double prev = guess, fPrev, fDerivPrev;
		double actual;
		int goodTries = 0;
		int badTries = 0;

		double previousDistance = 9999999.0;
		double distance;

		Polinom basePoli(*this);
		Polinom derivedPoli = GetDerived();

		int maxRoots = basePoli.GetDegree();

		do {
			fPrev = basePoli.GetValueAt(prev);
			fDerivPrev = derivedPoli.GetValueAt(prev);

			if (FloatHelper::IsZero(fDerivPrev)) {
				return result;
			}

			actual = prev - (fPrev / fDerivPrev);

			distance = fabs(actual - prev);

			if (distance < numeric_epsilon) {
				goodTries++;
			}
			else if (distance > previousDistance) {
				badTries++;
			}

			previousDistance = distance;
			prev = actual;

			if (badTries > max_wrong_tries) {
				return result;
			}

			if (goodTries > min_good_tries) {
				result.AddRoot(actual);
				basePoli = basePoli.divideByRoot(actual);
				derivedPoli = basePoli.GetDerived();
				goodTries = 0;
				badTries = 0;
				prev = guess;
				previousDistance = 9999999.0;
			}

		} while (result.GetNumberOfRoots() < maxRoots);

		return result;

	}

	RealRoots getRootsWithBisection(Polinom polinom, double initA, double initB, RealRoots *foundRoots) {
		RealRoots result;

		if (foundRoots != NULL) {
			result = *foundRoots;
		}

		int numOfSteps = 1000;

		double A = initA;
		double B = initB;
		double C;
		double step = (B - A) / numOfSteps;

		Polinom basePoli = polinom;

		int sigA = FloatHelper::Sig(basePoli.GetValueAt(A));
		int sigB = FloatHelper::Sig(basePoli.GetValueAt(B));
		int sigC;

		if (sigA == 0) {
			result.AddRoot(A);
			result = getRootsWithBisection(basePoli.divideByRoot(A), initA, initB, &result);
			return result;
		}

		if (sigB == 0) {
			result.AddRoot(B);
			result = getRootsWithBisection(basePoli.divideByRoot(B), initA, initB, &result);
			return result;
		}

		while (sigA == sigB && A < B) {
			A += step;
			sigA = FloatHelper::Sig(basePoli.GetValueAt(A));
		}

		if (!(A < B)) return result;

		// There is at least one root

		double root;
		do {
			if (sigA == 0) {
				root = A;
				break;
			}

			if (sigB == 0) {
				root = B;
				break;
			}

			C = 0.5 * (A + B);
			sigC = FloatHelper::Sig(basePoli.GetValueAt(C));

			if (sigC == 0) {
				root = C;
				break;
			}

			if (abs(A - B) < 0.00001) {
				root = C;
				break;
			}

			if (sigC == sigA) {
				A = C;
				sigA = sigC;
			}
			else if (sigC == sigB) {
				B = C;
				sigB = sigC;
			}

		} while (true);

		result.AddRoot(root);
		result = getRootsWithBisection(basePoli.divideByRoot(root), initA, initB, &result);
		return result;
	}


	friend ostream& operator<<(ostream& os, const Polinom& p);


};

const double Polinom::numeric_epsilon = 0.000001;
const int Polinom::max_wrong_tries = 10;
const int Polinom::min_good_tries = 3;

ostream& operator<<(ostream& os, const Polinom& p) {
	for (int i = p.degree; i >= 0; i--) {
		double next;
		os << (i == p.degree ? p.GetCoefficient(i) : fabs(p.GetCoefficient(i)));
		if (i > 0) {
			next = p.GetCoefficient(i - 1);
			os << " * x";
			if (i > 1) {
				os << "^" << i;
			}
			os << (next > 0.0 ? " + " : " - ");
		}
	}
	return os;
}

class Vector {
	double x;
	double y;
	double z;
	double h;

public:
	Vector(double x = 0.0, double y = 0.0, double z = 0.0, double h = 1.0) : x(x), y(y), z(z), h(h) {}

	double GetX() const { return x / h; }
	double GetY() const { return y / h; }
	double GetZ() const { return z / h; }
	double GetH() const { return h; }

	double GetRawX() const { return x; }
	double GetRawY() const { return y; }
	double GetRawZ() const { return z; }

	Vector& X(double newX) {
		x = newX;
		return *this;
	}
	Vector& Y(double newY) {
		y = newY;
		return *this;
	}
	Vector& Z(double newZ) {
		z = newZ;
		return *this;
	}
	Vector& H(double newH) {
		h = newH;
		return *this;
	}

	void Set(double vX = 0.0, double vY = 0.0, double vZ = 0.0, double vH = 1.0) {
		x = vX;
		y = vY;
		z = vZ;
		h = vH;
	}

	double GetLenght() { return (double)sqrt((x*x + y*y + z*z) / (h*h)); }

	Vector operator+(const Vector& v) { return Vector(GetX() + v.GetX(), GetY() + v.GetY(), GetZ() + v.GetZ(), 1.0); }
	Vector operator-(const Vector& v) { return Vector(GetX() - v.GetX(), GetY() - v.GetY(), GetZ() - v.GetZ(), 1.0); }
	Vector operator%(const Vector& v) { return Vector(GetY() * v.GetZ() - GetZ() * v.GetY(), GetZ() * v.GetX() - GetX() * v.GetZ(), GetX() * v.GetY() - GetY() * v.GetX()); }
	Vector operator*(const double& c) { return Vector(c * GetX(), c * GetY(), c * GetZ(), 1.0); }
	double operator*(const Vector& v) { return (GetX() * v.GetX() + GetY() * v.GetY() + GetZ() * v.GetZ()); }
	Vector operator/(const double& c) { return Vector(GetX() / c, GetY() / c, GetZ() / c, 1.0); }

	Vector& operator=(const Vector& v) {
		if (this != &v) {
			x = v.x;
			y = v.y;
			z = v.z;
			h = v.h;
		}
		return *this;
	}

	void Normalize() {
		double length = this->GetLenght();
		if (FloatHelper::IsZero(length)) {
			return;
		}

		double newX = (this->GetX()) / length;
		double newY = (this->GetY()) / length;
		double newZ = (this->GetZ()) / length;
		x = newX;
		y = newY;
		z = newZ;
		h = 1.0;
	}

	Vector GetNormalized() {
		Vector result(GetX(), GetY(), GetZ(), GetH());
		result.Normalize();
		return result;
	}
};


class Color {
	double r;
	double g;
	double b;

public:

	Color(double r = 0.0, double g = 0.0, double b = 0.0) : r(r), g(g), b(b) {}

	double GetR() const { return r; }
	double GetG() const { return g; }
	double GetB() const { return b; }

	Color& R(double newR) {
		r = newR;
		return *this;
	}

	Color& G(double newG) {
		g = newG;
		return *this;
	}

	Color& B(double newB) {
		b = newB;
		return *this;
	}

	void Set(double cR = 0.0, double cG = 0.0, double cB = 0.0) {
		r = cR;
		g = cG;
		b = cB;
	}

	Color& operator=(const Color& other) {
		if (this != &other) {
			r = other.r;
			g = other.g;
			b = other.b;
		}
		return *this;
	}

	Color operator+(const Color& other) { return Color(GetR() + other.GetR(), GetG() + other.GetG(), GetB() + other.GetB()); }
	Color operator-(const Color& other) { return Color(GetR() - other.GetR(), GetG() - other.GetG(), GetB() - other.GetB()); }
	Color operator*(const double& c) { return Color(GetR() * c, GetG() * c, GetB() * c); }
	Color operator/(const double& c) { return Color(GetR() / c, GetG() / c, GetB() / c); }
	Color operator*(const Color& other) { return Color(GetR() * other.GetR(), GetG() * other.GetG(), GetB() * other.GetB()); }
	Color operator/(const Color& other) { return Color(GetR() / other.GetR(), GetG() / other.GetG(), GetB() / other.GetB()); }

};

class Ray {
	Vector p0;
	Vector v;

public:
	Ray(Vector p0 = Vector(), Vector v = Vector()) : p0(p0), v(v) {}

	Ray& operator=(const Ray& other) {
		if (this != &other) {
			p0 = other.p0;
			v = other.v;
		}
		return *this;
	}

	Ray& P0(Vector newP0) {
		p0 = newP0;
		return *this;
	}

	Ray& V(Vector newV) {
		v = newV;
		return *this;
	}

	Ray Reflect(Vector point, Vector normal) {
		Vector baseDir = v.GetNormalized();
		Vector reflectedDir = baseDir - normal * (2.0 * (normal * baseDir));
		Ray reflected;
		reflected.P0(point + normal * EPSZILON_TOLAS).V(reflectedDir);
		return reflected;
	}

	Ray Refract(Vector point, Vector normal, double n) {
		Vector baseDir = v.GetNormalized();
		double cosalfa = normal * (baseDir * (-1.0));
		double sinalfanegyzet = 1.0 - cosalfa*cosalfa;
		double cosbetanegyzet = 1.0 - sinalfanegyzet / (n*n);
		double cosbeta = sqrt(cosbetanegyzet);
		Vector refractedDir = baseDir / n + normal * (cosalfa / n - cosbeta);
		Ray refracted;
		refracted.P0(point - normal * EPSZILON_TOLAS).V(refractedDir);
		return refracted;
	}

	Vector GetP0() { return p0; }
	Vector GetV() { return v; }

};

class LightSource {
	Color lightColor;
	double intensity;
	Vector lightPosition;
public:
	LightSource(
		double newIntensity = 0.0,
		Vector newPosition = Vector().X(0.0).Y(0.0).Z(0.0),
		Color newColor = Color().R(1.0).G(1.0).B(1.0)
	)
		:intensity(newIntensity), lightPosition(newPosition), lightColor(newColor) {}

	LightSource& LightColor(Color newColor) {
		lightColor = newColor;
		return *this;
	}

	Color GetLightColor() { return lightColor; }

	LightSource& Intensity(double newIntensity) {
		intensity = newIntensity;
		return *this;
	}

	double GetIntensity() { return intensity; }

	LightSource& LightPosition(Vector newPosition) {
		lightPosition = newPosition;
		return *this;
	}

	Vector getLightPosition() { return lightPosition; }

	Color GetLightAt(Vector point) {
		Vector distanceVector = lightPosition - point;
		double distance = distanceVector.GetLenght();
		double power;
		if (FloatHelper::IsZero(distance)) power = intensity;
		else power = intensity / (distance * distance);
		return lightColor * power;
	}


};

class Intersection {
	Vector intersectPoint;
	Vector normal;
	double t;
	int index;

public:
	Intersection() :t(-1.0), index(-1) {}

	Intersection& Param(double param) {
		t = param;
		return *this;
	}

	double GetParam() { return t; }

	Intersection& IntersectPoint(Vector mp) {
		intersectPoint = mp;
		return *this;
	}

	Vector GetIntersectPoint() { return intersectPoint; }

	Intersection& Normal(Vector norm) {
		normal = norm;
		return *this;
	}

	Vector GetNormal() { return normal; }

	Intersection& Index(int newIndex) {
		index = newIndex;
		return *this;
	}

	int GetIndex() { return index; }

	Intersection& operator=(const Intersection& other) {
		if (this != &other) {
			intersectPoint = other.intersectPoint;
			normal = other.normal;
			t = other.t;
			index = other.index;
		}
		return *this;
	}

};

class PhotonMap {
	Color cells[FOTONTERKEP_X][FOTONTERKEP_Y];
	bool postProcessDone;
public:

	PhotonMap() { postProcessDone = false; }

	Color GetPhoton(int x, int y) {
		if (!postProcessDone) PostProcess();
		return cells[x][y];
	}

	void AddPhoton(int x, int y, Color photonIntensity) {
		if ((x >= 0) && (x <= (FOTONTERKEP_X - 1)) && (y >= 0) && (y <= (FOTONTERKEP_Y - 1))) {
			cells[x][y] = cells[x][y] + photonIntensity;
			postProcessDone = false;
		}
	}

private:
	void PostProcess() {
		for (int a = 0; a < MAX_ATLAG; a++) {
			Color temp;
			for (int i = ATLAG_SUGAR; i < FOTONTERKEP_X - ATLAG_SUGAR; i++) {
				for (int j = ATLAG_SUGAR; j < FOTONTERKEP_Y - ATLAG_SUGAR; j++) {
					temp.Set(0.0, 0.0, 0.0);
					for (int x = i - ATLAG_SUGAR; x <= i + ATLAG_SUGAR; x++) {
						for (int y = j - ATLAG_SUGAR; y <= j + ATLAG_SUGAR; y++) {
							Color original = cells[x][y];
							temp = temp + cells[x][y];
						}
					}

					temp = temp / ((2 * ATLAG_SUGAR + 1) * (2 * ATLAG_SUGAR + 1));
					cells[i][j] = temp;
				}
			}
		}
		postProcessDone = true;
	}

public:

	void ReadFromFile(const char *fileName) {
		// Fájlból feltölteni a fotontérképet

		ifstream file(fileName);

		if (!file.is_open()) return;

		double R, G, B;

		for (int i = 0; i < FOTONTERKEP_X; i++) {
			for (int j = 0; j < FOTONTERKEP_Y; j++) {
				file >> R >> G >> B;
				cells[i][j].R(R).G(G).B(B);
			}
		}

		file >> postProcessDone;

		file.close();
	}

	void WriteToFile(const char *fileName) {
		// Fotontérkép fájlbaírása

		ofstream file(fileName);

		if (!file.is_open()) return;

		for (int i = 0; i < FOTONTERKEP_X; i++) {
			for (int j = 0; j < FOTONTERKEP_Y; j++) {
				file << cells[i][j].GetR() << cells[i][j].GetG() << cells[i][j].GetB();
			}
		}

		file << postProcessDone;

		file.close();
	}

};

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
	Color kd, ks, ka, F0, n, k;
	double shine;
	ProceduralMode proceduralMode;
	PhotonMap* photonMap;

public:

	Object() {
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

	Object& Diffuse(bool diffuse = true) {
		isDiffuse = diffuse;
		return *this;
	}

	bool IsDiffuse() { return isDiffuse; }

	Object& Reflective(bool reflective = true) {
		isReflective = reflective;
		return *this;
	}

	bool IsReflective() { return isReflective; }

	Object& Refractive(bool refractive = true) {
		isRefractive = refractive;
		return *this;
	}

	bool IsRefractive() { return isRefractive; }

	Object& Kd(Color newKd) {
		kd = newKd;
		return *this;
	}

	Color GetKd() { return kd; }

	Object& Ks(Color newKs) {
		ks = newKs;
		return *this;
	}

	Color GetKs() { return ks; }

	Object& Ka(Color newKa) {
		ka = newKa;
		return *this;
	}

	Color GetKa() { return ka; }

	Object& Shine(double newShine) {
		shine = newShine;
		return *this;
	}

	double GetShine() { return shine; }

	void calculateFresnel() {
		F0
			.R(((n.GetR() - 1.0)*(n.GetR() - 1.0) + k.GetR()*k.GetR()) / ((n.GetR() + 1.0)*(n.GetR() + 1.0) + k.GetR()*k.GetR()))
			.G(((n.GetG() - 1.0)*(n.GetG() - 1.0) + k.GetG()*k.GetG()) / ((n.GetG() + 1.0)*(n.GetG() + 1.0) + k.GetG()*k.GetG()))
			.B(((n.GetB() - 1.0)*(n.GetB() - 1.0) + k.GetB()*k.GetB()) / ((n.GetB() + 1.0)*(n.GetB() + 1.0) + k.GetB()*k.GetB()));
	}

	Object& N(Color newN) {
		n = newN;
		calculateFresnel();
		return *this;
	}

	Color GetN() { return n; }

	Object& K(Color newK) {
		k = newK;
		calculateFresnel();
		return *this;
	}

	Color GetK() { return k; }

	Color GetFresnel(Vector N, Vector V) {
		double cosalfa = fabs(N*V);
		Color result;
		Color a;
		Color b;
		double factor = pow((1.0 - cosalfa), 5.0);
		a = F0;
		b.R(1.0 - F0.GetR()).G(1.0 - F0.GetG()).B(1.0 - F0.GetB());
		result = a + b*factor;
		return result;
	}

	Object& ProcMode(ProceduralMode newMode) {
		proceduralMode = newMode;
		return *this;
	}

	ProceduralMode GetProcMode() { return proceduralMode; }

	void SetPhotonMap(PhotonMap *newMap) {
		photonMap = newMap;
	}

	PhotonMap* GetPhotonMap() {
		return photonMap;
	}

	void WritePhotonMapToFile(const char* fileName) {
		if (photonMap != NULL) photonMap->WriteToFile(fileName);
	}

	void ReadPhotonMapFromFile(const char* fileName) {
		if (photonMap != NULL) photonMap->ReadFromFile(fileName);
	}

	virtual Intersection Intersect(Ray ray) = 0;

	virtual void AddPhoton(Color intensity, Vector position) {};

	virtual Color GetPhoton(Vector position) {
		return Color().R(0.0).G(0.0).B(0.0);
	}

	virtual Color GetProceduralColor(Vector position) {
		return Color().R(1.0).G(1.0).B(1.0);
	}

};

class Cylinder :public Object {
	Vector basePoint;
	double radius;
	double height;
	bool solid;
public:
	Cylinder() :Object() { solid = false; }

	Cylinder& BasePoint(Vector newBasePoint) {
		basePoint = newBasePoint;
		return *this;
	}

	Vector GetBasePoint() { return basePoint; }

	Cylinder& Radius(double newRadius) {
		radius = newRadius;
		return *this;
	}

	double GetRadius() { return radius; }

	Cylinder& Height(double newHeight) {
		height = newHeight;
		return *this;
	}

	double GetHeight() { return height; }

	Cylinder& Solid(bool newSolid = true) {
		solid = newSolid;
		return *this;
	}

	bool IsSolid() { return solid; }

	Intersection Intersect(Ray ray) {
		Intersection result;
		result.Param(-1.0);

		double px = ray.GetP0().GetX();
		double py = ray.GetP0().GetY();
		double vx = ray.GetV().GetX();
		double vy = ray.GetV().GetY();
		double hx = basePoint.GetX();
		double hy = basePoint.GetY();
		double r = radius;

		double a = vx*vx + vy*vy;
		double b = 2.0*(px*vx + py*vy - vx*hx - vy*hy);
		double c = px*px + py*py + hx*hx + hy*hy - 2.0*(px*hx + py*hy) - r*r;
		double discr = b*b - 4 * a*c;
		if (discr < 0.0) return result;

		double bottomHeight = basePoint.GetZ();
		double topHeight = bottomHeight + height;
		if (FloatHelper::IsZero(discr)) {
			double t = (-1.0*b) / (2.0*a);
			Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
			Vector helper;
			double z = intersectPoint.GetZ();
			if ((z < bottomHeight) ||
				(z > topHeight)
				|| (t < 0.0)) {
				return result;
			}
			else {
				helper = intersectPoint;
				helper.Z(bottomHeight);
				Vector normal = (helper - basePoint).GetNormalized();
				result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
				return result;
			}
		}
		else {
			double t1 = (-1.0*b - sqrt(discr)) / (2.0*a);
			double t2 = (-1.0*b + sqrt(discr)) / (2.0*a);

			Vector intersectPoint1 = ray.GetP0() + ray.GetV()*t1;
			Vector intersectPoint2 = ray.GetP0() + ray.GetV()*t2;
			Vector helper;
			double z1 = intersectPoint1.GetZ();
			double z2 = intersectPoint2.GetZ();

			if ((z1 >= bottomHeight) && (z1 <= topHeight) && (t1 > 0.0)) {
				// Elsõ metszéspont magassága stimmel
				helper = intersectPoint1;
				helper.Z(bottomHeight);
				Vector normal = (helper - basePoint).GetNormalized();
				result.Param(t1).IntersectPoint(intersectPoint1).Normal(normal);
				return result;
			}

			if (solid) {
				// Tömör hengernél ellenõrizni kell a tetõ és alaplap metszését
				double pz = ray.GetP0().GetZ();
				double vz = ray.GetV().GetZ();

				int ABOVE = 1;
				int CORRECT = 2;
				int BELOW = 3;

				int firstHeight, secondHeight;

				if (z1 > topHeight) {
					firstHeight = ABOVE;
				}
				else if (z1 < bottomHeight) {
					firstHeight = BELOW;
				}
				else {
					firstHeight = CORRECT;
				}

				if (z2 > topHeight) {
					secondHeight = ABOVE;
				}
				else if (z2 < bottomHeight) {
					secondHeight = BELOW;
				}
				else {
					secondHeight = CORRECT;
				}

				if (firstHeight == ABOVE && (secondHeight == CORRECT || secondHeight == BELOW)) {
					// Tetõlap metszve

					// pz + vz*t = teto => t = (teto - pz) / vz
					double t = (topHeight - pz) / vz;
					Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
					Vector normal(0.0, 0.0, 1.0);
					result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
					return result;
				}

				if (firstHeight == BELOW && (secondHeight == CORRECT || secondHeight == ABOVE)) {
					// Alaplap metszve

					// pz + vz*t = alapColort => t = (alapColort - pz) / vz
					double t = (bottomHeight - pz) / vz;
					Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
					Vector normal(0.0, 0.0, -1.0);
					result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
					return result;
				}

			}

			// Nem tömör hengernél eddigi módszer
			if ((z2 > bottomHeight) && (z2 < topHeight) && (t2 > 0.0)) {
				helper = intersectPoint2;
				helper.Z(bottomHeight);
				Vector normal = (helper - basePoint).GetNormalized();
				normal = normal * -1.0;
				result.Param(t2).IntersectPoint(intersectPoint2).Normal(normal);
				return result;
			}

			return result;
		}
	}
};

class Ellipsoid :public Object {
	Vector center;
	double _A, _B, _C;
public:
	Ellipsoid() :Object() {}

	Ellipsoid& Center(Vector newCenter) {
		center = newCenter;
		return *this;
	}

	Vector GetCenter() { return center; }

	Ellipsoid& A(double newA) {
		_A = newA;
		return *this;
	}

	double GetA() { return _A; }

	Ellipsoid& B(double newB) {
		_B = newB;
		return *this;
	}

	double GetB() { return _B; }

	Ellipsoid& C(double newC) {
		_C = newC;
		return *this;
	}

	double GetC() { return _C; }

	Intersection Intersect(Ray ray) {
		Intersection result;

		result.Param(-1.0);

		Vector i(1.0, 0.0, 0.0);
		Vector j(0.0, 1.0, 0.0);
		Vector k(0.0, 0.0, 1.0);

		double px = ray.GetP0().GetX();
		double py = ray.GetP0().GetY();
		double pz = ray.GetP0().GetZ();
		double vx = ray.GetV().GetX();
		double vy = ray.GetV().GetY();
		double vz = ray.GetV().GetZ();
		double ex = center.GetX();
		double ey = center.GetY();
		double ez = center.GetZ();

		double eh1 = _B*_B * _C*_C;
		double eh2 = _A*_A * _C*_C;
		double eh3 = _A*_A * _B*_B;
		double R = _A*_A *_B*_B* _C*_C;

		double a = eh1*vx*vx + eh2*vy*vy + eh3*vz*vz;
		double b = 2.0 * (eh1 * (px*vx - vx*ex) + eh2 * (py*vy - vy*ey) + eh3 * (pz*vz - vz*ez));
		double c = eh1 * (px*px + ex*ex - 2.0 * px*ex) + eh2 * (py*py + ey*ey - 2.0 * py*ey) + eh3 * (pz*pz + ez*ez - 2.0 * pz*ez) - R;

		double discr = b*b - 4 * a*c;

		if (discr < 0.0) return result;

		double t1 = (-1.0*b - sqrt(discr)) / (2.0*a);
		double t2 = (-1.0*b + sqrt(discr)) / (2.0*a);

		if ((t1 > 0.0) || (t2 > 0.0)) {
			double t = (t1 > 0.0) ? t1 : t2;
			Vector intersectPoint = ray.GetP0() + ray.GetV()*t;
			double ix = intersectPoint.GetX();
			double iy = intersectPoint.GetY();
			double iz = intersectPoint.GetZ();
			double nx = 2.0*(ix - ex) / (_A*_A);
			double ny = 2.0*(iy - ey) / (_B*_B);
			double nz = 2.0*(iz - ez) / (_C*_C);
			Vector normal = (i*nx + j*ny + k*nz).GetNormalized();
			result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
		}

		return result;
	}
};

class Rect :public Object {
	Vector corner1;
	Vector corner2;
	Vector corner3;
	Vector corner4;
public:
	Rect() : Object() {}

	Rect& Corner1(Vector newCorner1) {
		corner1 = newCorner1;
		return *this;
	}

	Vector GetCorner1() { return corner1; }

	Rect& Corner2(Vector newCorner2) {
		corner2 = newCorner2;
		return *this;
	}

	Vector GetCorner2() { return corner2; }

	Rect& Corner3(Vector newCorner3) {
		corner3 = newCorner3;
		return *this;
	}

	Vector GetCorner3() { return corner3; }

	Rect& Corner4(Vector newCorner4) {
		corner4 = newCorner4;
		return *this;
	}

	Vector GetCorner4() { return corner4; }

	double GetSizeX() { return (corner2 - corner1).GetLenght(); }
	double GetSizeY() { return (corner3 - corner1).GetLenght(); }

	Intersection Intersect(Ray ray) {
		Intersection result;

		Vector normal = ((corner2 - corner1) % (corner3 - corner1)).GetNormalized();

		Vector p = ray.GetP0();
		Vector v = ray.GetV();

		double t = (corner1*normal - p*normal) / (v*normal);

		if (t >= 0.0) {
			Vector intersectPoint = p + v*t;

			double test1 = ((corner2 - corner1) % (intersectPoint - corner1)) * normal;
			double test2 = ((corner3 - corner2) % (intersectPoint - corner2)) * normal;
			double test3 = ((corner4 - corner3) % (intersectPoint - corner3)) * normal;
			double test4 = ((corner1 - corner4) % (intersectPoint - corner4)) * normal;

			bool good1 = (test1 > 0.0);
			bool good2 = (test2 > 0.0);
			bool good3 = (test3 > 0.0);
			bool good4 = (test4 > 0.0);

			if (good1 && good2 && good3 && good4) {
				result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
			}
		}

		return result;
	}

	virtual void AddPhoton(Color intensity, Vector position) {
		if (photonMap != NULL) {
			Vector a = (corner2 - corner1).GetNormalized();
			Vector c = position - corner1;
			double e = c.GetLenght();

			double photonCellSizeX = GetSizeX() / FOTONTERKEP_X;
			double photonCellSizeY = GetSizeY() / FOTONTERKEP_Y;

			int photonX = 0;
			int	photonY = 0;

			if (e > 0.0) {
				c.Normalize();
				double cosalfa = c*a;
				double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

				double pointX = e*cosalfa;
				double pointY = e*sinalfa;

				while (pointX >= photonCellSizeX) {
					photonX++;
					pointX -= photonCellSizeX;
				}

				while (pointY >= photonCellSizeY) {
					photonY++;
					pointY -= photonCellSizeY;
				}
			}

			photonMap->AddPhoton(photonX, photonY, intensity);
		}

	}

	virtual Color GetPhoton(Vector position) {
		Color photon = Object::GetPhoton(position);

		if (photonMap != NULL) {
			Vector a = (corner2 - corner1).GetNormalized();
			Vector c = position - corner1;
			double e = c.GetLenght();

			double photonCellSizeX = GetSizeX() / FOTONTERKEP_X;
			double photonCellSizeY = GetSizeY() / FOTONTERKEP_Y;

			int photonX = 0;
			int	photonY = 0;

			if (e > 0.0) {
				c.Normalize();
				double cosalfa = c*a;
				double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

				double pointX = e*cosalfa;
				double pointY = e*sinalfa;

				while (pointX >= photonCellSizeX) {
					photonX++;
					pointX -= photonCellSizeX;
				}

				while (pointY >= photonCellSizeY) {
					photonY++;
					pointY -= photonCellSizeY;
				}
			}

			photon = photonMap->GetPhoton(photonX, photonY);
		}

		return photon;
	}

	virtual Color GetProceduralColor(Vector position) {
		Color result = Object::GetProceduralColor(position);

		if (TEXTURA) {
			int i = 0;
			int j = 0;

			Vector a = (corner2 - corner1).GetNormalized();

			double unitX = GetSizeX() / TEXTURA_X;
			double unitY = GetSizeY() / TEXTURA_Y;

			Vector c = position - corner1;
			double e = c.GetLenght();

			if (e > 0.0) {
				c.Normalize();
				double cosalfa = c*a;
				double sinalfa = sqrt(1.0 - cosalfa*cosalfa);

				double diffX = e*cosalfa;
				double diffY = e*sinalfa;

				while (diffX > unitX) {
					diffX -= unitX;
					i++;
				}
				while (diffY > unitY) {
					diffY -= unitY;
					j++;
				}
			}

			double factor = sqrt((double)(i^j)) + sqrt((double)(j^i));
			int multiplier = (int)(floor(factor)) % 5;

			result = result * (0.1 + 0.1*multiplier);
		}

		return result;
	}

};

class Cube : public Object {
	Rect front;
	Rect back;
	Rect left;
	Rect right;
	Rect top;
	Rect bottom;

public:
	Cube() {
		Vector FLD(-2.0, -2.0, -2.0);
		Vector FRD(2.0, -2.0, -2.0);
		Vector FRU(2.0, -2.0, 2.0);
		Vector FLU(-2.0, -2.0, 2.0);
		Vector BLD(-2.0, 2.0, -2.0);
		Vector BRD(2.0, 2.0, -2.0);
		Vector BRU(2.0, 2.0, 2.0);
		Vector BLU(-2.0, 2.0, 2.0);

		front.Corner1(FLD).Corner2(FRD).Corner3(FRU).Corner4(FLU);
		back.Corner1(BRD).Corner2(BLD).Corner3(BLU).Corner4(BRU);
		left.Corner1(BLD).Corner2(FLD).Corner3(FLU).Corner4(BLU);
		right.Corner1(FRD).Corner2(BRD).Corner3(BRU).Corner4(FRU);
		top.Corner1(FLU).Corner2(FRU).Corner3(BRU).Corner4(BLU);
		bottom.Corner1(FLD).Corner2(BLD).Corner3(BRD).Corner4(FRD);
	}

	Intersection Intersect(Ray ray) {
		Intersection result;
		Intersection actual;

		double tmin = 99999999999.0;
		actual = front.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		actual = back.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		actual = left.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		actual = right.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		actual = top.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		actual = bottom.Intersect(ray);
		if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
			tmin = actual.GetParam();
			result = actual;
		}
		return result;
	}

};

class Triangle :public Object {
	Vector vertex1;
	Vector vertex2;
	Vector vertex3;

public:
	Triangle() : Object() {}

	Triangle& Vertex1(Vector newVertex1) {
		vertex1 = newVertex1;
		return *this;
	}

	Vector GetVertex1() { return vertex1; }

	Triangle& Vertex2(Vector newVertex2) {
		vertex2 = newVertex2;
		return *this;
	}

	Vector GetVertex2() { return vertex2; }

	Triangle& Vertex3(Vector newVertex3) {
		vertex3 = newVertex3;
		return *this;
	}

	Vector GetVertex3() { return vertex3; }


	Intersection Intersect(Ray ray) {
		Intersection result;

		Vector normal = ((vertex2 - vertex1) % (vertex3 - vertex1)).GetNormalized();

		Vector p = ray.GetP0();
		Vector v = ray.GetV();

		double t = (vertex1*normal - p*normal) / (v*normal);

		if (t > 0.0) {
			Vector intersectPoint = p + v*t;

			double test1 = ((vertex2 - vertex1) % (intersectPoint - vertex1)) * normal;
			double test2 = ((vertex3 - vertex2) % (intersectPoint - vertex2)) * normal;
			double test3 = ((vertex1 - vertex3) % (intersectPoint - vertex3)) * normal;

			bool good1 = (test1 > 0.0);
			bool good2 = (test2 > 0.0);
			bool good3 = (test3 > 0.0);

			if (good1 && good2 && good3) {
				result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
			}
		}

		return result;
	}

	void Beallit(Vector cs1, Vector cs2, Vector cs3) {
		vertex1 = cs1;
		vertex2 = cs2;
		vertex3 = cs3;
	}

};

class TriangleMesh : public Object {
	Triangle* triangles[MAX_HAROMSZOG_SZAM];
	int numOfTriangles;
public:
	TriangleMesh() : Object() {
		Clear();
	}

	int GetNumberOfTriangles() { return numOfTriangles; }

	void Clear() {
		numOfTriangles = 0;
		for (int i = 0; i<MAX_HAROMSZOG_SZAM; i++) triangles[i] = NULL;
	}

	void AddTriangle(Triangle *newTriangle) {
		if (numOfTriangles < MAX_HAROMSZOG_SZAM && newTriangle != NULL) {
			triangles[numOfTriangles++] = newTriangle;
		}
	}

	Intersection Intersect(Ray ray) {
		Intersection result;
		Intersection actual;

		double tmin = 99999999999.0;
		if (numOfTriangles > 0) {
			for (int i = 0; i < numOfTriangles; i++) {
				actual = triangles[i]->Intersect(ray);
				if ((actual.GetParam() > 0.0) && (actual.GetParam() < tmin)) {
					tmin = actual.GetParam();
					result = actual;
					result.Index(i);
				}
			}
		}
		return result;
	}

};

class Diamond :public Object {
	Ellipsoid *boundingVolume;
	TriangleMesh *diamondBody;
public:

	Diamond() {
		boundingVolume = NULL;
		diamondBody = NULL;
	}

	Diamond& BoundingVolume(Ellipsoid *newBounding) {
		boundingVolume = newBounding;
		return *this;
	}

	Ellipsoid* GetBoundingVolume() { return boundingVolume; }

	Diamond& DiamondBody(TriangleMesh *newBody) {
		diamondBody = newBody;
		return *this;
	}

	TriangleMesh* GetDiamondBody() { return diamondBody; }

	Intersection Intersect(Ray ray) {
		Intersection result;

		if (boundingVolume != NULL && diamondBody != NULL) {
			Intersection boundingIntersection = boundingVolume->Intersect(ray);

			if (boundingIntersection.GetParam() >= 0.0) {
				result = diamondBody->Intersect(ray);
			}
		}

		return result;
	}

};

class HeartShape : public Object {

	Cube cube;

	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;

public:
	HeartShape() {
		minX = 99999999.0;
		maxX = -99999999.0;
		minY = 99999999.0;
		maxY = -99999999.0;
		minZ = 99999999.0;
		maxZ = -99999999.0;
	}

	Intersection Intersect(Ray ray) {
		Intersection result;

		Intersection test;
		test = cube.Intersect(ray);
		double tmin = test.GetParam();

		if (tmin < 0.0) return result;

		Vector newP0 = test.GetIntersectPoint() + ray.GetV() * 0.001;
		Ray ray2(newP0, ray.GetV());

		Intersection test2 = cube.Intersect(ray2);
		double ttemp = test2.GetParam();

		if (ttemp < 0.0) return result;

		Vector testIntersect = test2.GetIntersectPoint();

		double tmax = -1.0;

		double tix = testIntersect.GetX();
		double tiy = testIntersect.GetY();
		double tiz = testIntersect.GetZ();

		double px = ray.GetP0().GetX();
		double py = ray.GetP0().GetY();
		double pz = ray.GetP0().GetZ();

		double vx = ray.GetV().GetX();
		double vy = ray.GetV().GetY();
		double vz = ray.GetV().GetZ();

		if (!FloatHelper::IsZero(vx)) {
			tmax = (tix - px) / vx;
		}
		else if (!FloatHelper::IsZero(vy)) {
			tmax = (tiy - py) / vy;
		}
		else if (!FloatHelper::IsZero(vy)) {
			tmax = (tiz - pz) / vz;
		}

		if (tmax < 0.0) return result;

		//double guess = test.GetParam();

		double c1 = vx * vx + 2.25 * vy * vy + vz * vz;
		double c2 = 2 * px * vx + 4.5 * py * vy + 2 * pz * vz;
		double c3 = px * px + 2.25 * py * py + pz * pz - 1;

		double c4 = c1 * c1;
		double c5 = 2 * c1 * c2;
		double c6 = c2 * c2 + 2 * c1 * c3;
		double c7 = 2 * c2 * c3;
		double c8 = c3 * c3;

		double c9 = c1 * c4;
		double c10 = c1 * c5 + c2 * c4;
		double c11 = c1 * c6 + c2 * c5 + c3 * c4;
		double c12 = c1 * c7 + c2 * c6 + c3 * c5;
		double c13 = c1 * c8 + c2 * c7 + c3 * c6;
		double c14 = c2 * c8 + c3 * c7;
		double c15 = c3 * c8;

		double c16 = vy * vy;
		double c17 = 2 * py * vy;
		double c18 = py * py;

		double c19 = vz * vz * vz;
		double c20 = 3 * pz * vz * vz;
		double c21 = 3 * pz * pz * vz;
		double c22 = pz * pz * pz;

		double c23 = c16 * c19;
		double c24 = c16 * c20 + c17 * c19;
		double c25 = c16 * c21 + c17 * c20 + c18 * c19;
		double c26 = c16 * c22 + c17 * c21 + c18 * c20;
		double c27 = c17 * c22 + c18 * c21;
		double c28 = c18 * c22;

		double c29 = vx * vx;
		double c30 = 2 * px * vx;
		double c31 = px * px;

		double c32 = c29 * c19;
		double c33 = c29 * c20 + c30 * c19;
		double c34 = c29 * c21 + c30 * c20 + c31 * c19;
		double c35 = c29 * c22 + c30 * c21 + c31 * c20;
		double c36 = c30 * c22 + c31 * c21;
		double c37 = c31 * c22;

		double constant1 = -9.0 / 80;

		double constant2 = -1.0;

		double a6 = c9;
		double a5 = c10 + constant1 * c23 + constant2 * c32;
		double a4 = c11 + constant1 * c24 + constant2 * c33;
		double a3 = c12 + constant1 * c25 + constant2 * c34;
		double a2 = c13 + constant1 * c26 + constant2 * c35;
		double a1 = c14 + constant1 * c27 + constant2 * c36;
		double a0 = c15 + constant1 * c28 + constant2 * c37;

		Polinom heartPoli(6);

		heartPoli.SetCoefficient(6, a6);
		heartPoli.SetCoefficient(5, a5);
		heartPoli.SetCoefficient(4, a4);
		heartPoli.SetCoefficient(3, a3);
		heartPoli.SetCoefficient(2, a2);
		heartPoli.SetCoefficient(1, a1);
		heartPoli.SetCoefficient(0, a0);

		//RealRoots roots = heartPoli.getRootsWithNewton(10.0);
		RealRoots roots = heartPoli.getRootsWithBisection(heartPoli, tmin, tmax, NULL);

		/*
		ofstream newtonFile("newton.txt", ios::out | ios::app);
		if ( newtonFile.is_open() ) {
		newtonFile << "---------------------------------------" << endl;
		newtonFile << "---------------------------------------" << endl;
		newtonFile << heartPoli << " solutions: " << endl;
		newtonFile << "---------------------------------------" << endl;
		newtonFile << roots << endl;
		newtonFile << "---------------------------------------" << endl;
		newtonFile << "---------------------------------------" << endl;
		newtonFile.close();
		}
		*/

		int number = roots.GetNumberOfRoots();

		if (number == 0) {
			return result;
		}

		double t = roots.GetMinPositiveRootIfAny();

		if (t <= 0.0) {
			return result;
		}

		Vector intersectPoint = ray.GetP0() + ray.GetV() * t;

		double ix = intersectPoint.GetX();
		double iy = intersectPoint.GetY();
		double iz = intersectPoint.GetZ();

		if (ix < minX) minX = ix;
		if (ix > maxX) maxX = ix;
		if (iy < minY) minY = iy;
		if (iy > maxY) maxY = iy;
		if (iz < minZ) minZ = iz;
		if (iz > maxZ) maxZ = iz;

		double n1 = ix * ix + 2.25f * iy * iy + iz * iz - 1.0;
		double n2 = n1 * n1;

		double parcX = 6 * ix * n2 - 2 * ix * iz * iz * iz;
		double parcY = 13.5 * iy * n2 - (9.0 / 40) * iy * iz * iz * iz;
		double parcZ = 6 * iz * n2 - 3 * ix * ix * iz * iz - (27.0 / 80) * iy * iy * iz * iz;

		/*
		double dxdy = parcX / parcY;
		double dxdz = parcX / parcZ;
		double dydx = parcY / parcX;
		double dydz = parcY / parcZ;
		double dzdx = parcZ / parcX;
		double dzdy = parcZ / parcY;

		double nx = 2 * n2 * (2 * ix + 2.25 * 2 * iy * dydx + 2 * iz * dzdx ) - (2 * ix * iz * iz * iz + ix * ix * 3 * iz * iz * dzdx ) - (9.0/80) * (2*iy*dydx * iz*iz*iz + iy*iy*3*iz*iz*dzdx );
		double ny = 2 * n2 * (2 * ix * dxdy + 2.25 * 2 * iy + 2 * iz * dzdy) - (2 * ix * dxdy * iz * iz * iz + ix * ix * 3 * iz * iz * dzdy ) - ( 9.0/80 ) * (2* iy * iz * iz * iz + iy * iy * 3 * iz * iz * dzdy );
		double nz = 2 * n2 * (2 * ix * dxdz + 2.25 * 2 * iy * dydz + 2 * iz ) - (2 * ix * dxdz * iz * iz * iz + ix * ix * 3 * iz * iz ) - (9.0 / 80 ) * (2 * iy * dydz * iz * iz * iz + iy * iy * 3 * iz * iz);

		Vector normal(nx, ny, nz);
		*/
		Vector normal(parcX, parcY, parcZ);
		normal.Normalize();
		result.Param(t).IntersectPoint(intersectPoint).Normal(normal);
		return result;
	}
};

class Camera {
	Vector lookAtPoint;
	Vector right;
	Vector up;
	Vector eye;
	double fovDegree;

public:
	Camera() {
		lookAtPoint.X(0.0).Y(0.0).Z(0.0).H(1.0);
		right.X(0.0).Y(0.0).Z(0.0).H(1.0);
		up.X(0.0).Y(0.0).Z(0.0).H(1.0);
		eye.X(0.0).Y(0.0).Z(0.0).H(1.0);
		fovDegree = 0.0;
	}

	Camera& LookAt(Vector newLookAt) {
		lookAtPoint = newLookAt;
		CalculateEye();
		return *this;
	}

	Vector GetLookAtPoint() { return lookAtPoint; }

	Camera& Right(Vector newRight) {
		right = newRight;
		CalculateEye();
		return *this;
	}

	Vector GetRightVector() { return right; }

	Camera& Up(Vector newUp) {
		up = newUp;
		CalculateEye();
		return *this;
	}

	Vector GetupVector() { return up; }

	Camera& FovDegree(double newFov) {
		fovDegree = newFov;
		CalculateEye();
		return *this;
	}

	double GetFovDegree() { return fovDegree; }

	Vector GetEye() { return eye; }

private:
	void CalculateEye() {
		double halfFovRad = ((fovDegree / 2.0) / 180.0) * PI;
		double tanHalfFov = tan(halfFovRad);
		double eyeDistance = (double)KEPMAGASSAG / (2.0 * tanHalfFov);
		Vector cameraAxisNormal = (right % up).GetNormalized();
		eye = lookAtPoint + cameraAxisNormal * eyeDistance;
	}

	Vector CenterOfPixel(int x, int z) {
		Vector result;
		//result = lookAtPoint + right * ( ((double)2*x/KEPSZELESSEG) - 1.0 ) + up * ( ((double)2*z/KEPMAGASSAG)-1 );
		result = lookAtPoint + right * (((double)2 * x / KEPSZELESSEG) - 1.0) + up * (1.0 - ((double)(2 * z - 2) / KEPMAGASSAG));
		return result;
	}

public:
	Ray GetRay(int x, int z) {
		Ray result;
		result.P0(eye).V(CenterOfPixel(x, z) - eye);
		return result;
	}

};

class Scene {
	Color ambient;
	Camera* camera;
	Object* objects[MAX_OBJ_SZAM];
	int numOfObjects;
	LightSource* lights[MAX_FENY_SZAM];
	int numOfLights;

public:
	Scene() {
		ambient.R(0.0).G(0.0).B(0.0);
		camera = NULL;
		ClearObjects();
		ClearLights();
	}

	Scene& Ambient(Color newAmbient) {
		ambient = newAmbient;
		return *this;
	}

	Color GetAmbient() { return ambient; }

	Scene& Cam(Camera *newCamera) {
		camera = newCamera;
		return *this;
	}

	Camera* GetCam() { return camera; }

	Scene& AddObject(Object *newObject) {
		if (numOfObjects < MAX_OBJ_SZAM) {
			objects[numOfObjects++] = newObject;
		}
		return *this;
	}

	int GetNumberOfObjects() { return numOfObjects; }

	void ClearObjects() {
		for (int i = 0; i < MAX_OBJ_SZAM; i++) objects[i] = NULL;
		numOfObjects = 0;
	}

	Scene& AddLightSource(LightSource *newLight) {
		if (numOfLights < MAX_FENY_SZAM) {
			lights[numOfLights++] = newLight;
		}
		return *this;
	}

	int GetNumberOfLightSources() { return numOfLights; }

	void ClearLights() {
		for (int i = 0; i < MAX_FENY_SZAM; i++) lights[i] = NULL;
		numOfLights = 0;
	}

	Intersection GetFirstIntersect(Ray ray) {
		Intersection result;
		Intersection actual;
		double tmin = 99999999999.0;
		if (numOfObjects > 0) {
			Intersection temp;
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

	Color GetDirectLight(Vector point, Vector N, Vector V, int index) {
		// point = point on surface
		// N = normal on surface
		// V = view direction = vector from point to the eye of the camera
		// index = index of the object thats lighened color we are calculating

		// outerPoint = point + epsilon * normal
		Vector outerPoint = point + N * EPSZILON_TOLAS;

		Object *object = objects[index];

		Color ka = object->GetKa();
		Color kd = object->GetKd();
		Color ks = object->GetKs();
		double shine = object->GetShine();
		Color procedural = object->GetProceduralColor(point);
		Object::ProceduralMode proceduralMode = object->GetProcMode();

		switch (proceduralMode)
		{
		case Object::ADD:
			kd = kd + procedural;
			ka = ka + procedural * PI;
			break;
		case Object::MULTIPLY:
			kd = kd * procedural;
			ka = ka * procedural * PI;
			break;
		case Object::REPLACE:
			kd = procedural;
			ka = procedural * PI;
			break;
		case Object::NONE:
		default:
			break;
		}

		// result = Ka * La
		Color result = ka * ambient;

		Vector lightPosition;
		// L = direction of shadow ray = vector from outerPoint to actual lightsource
		Vector L;
		double lightDistance;
		Intersection lightIntersection;
		Vector intersectVector;

		LightSource *actualLightSource;

		for (int i = 0; i < numOfLights; i++) {
			actualLightSource = lights[i];
			lightPosition = actualLightSource->getLightPosition();
			L = lightPosition - outerPoint;
			lightDistance = L.GetLenght();

			Ray shadowRay(outerPoint, L);
			lightIntersection = GetFirstIntersect(shadowRay);
			intersectVector = lightIntersection.GetIntersectPoint() - point;

			if ((lightIntersection.GetIndex() == -1) || (intersectVector.GetLenght() > lightDistance)) {
				L.Normalize();
				Vector H = (L + V);
				H.Normalize();
				double costheta = N * L;
				double cosdelta = N * H;

				// Lin: Light Source's color in point
				Color Lin = actualLightSource->GetLightAt(point);

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

	Color Trace(Ray ray, int depth, bool inside) {
		if (depth > D_MAX) return ambient;

		Intersection firstIntersection = GetFirstIntersect(ray);
		int index = firstIntersection.GetIndex();
		if (index == -1) return ambient;

		Object *object = objects[index];

		Vector N = firstIntersection.GetNormal();
		// IntersectPoint = x
		Vector point = firstIntersection.GetIntersectPoint();
		Vector V = (camera->GetEye() - point).GetNormalized();
		// RAy DIR = v
		Vector rayDir = ray.GetV().GetNormalized();
		if (N * rayDir > 0) N = N * -1.0;

		// Result = 0
		Color result;

		// If diffuse: Result = Result + DirectLight
		if (object->IsDiffuse() && LAMPAFENY) {
			bool LOL;
			if (depth > 0) {
				LOL = true;
			}
			result = result + GetDirectLight(point, N, V, index);
		}

		Color objectN;
		Color fresnel;
		Color kt;
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

	void Shoot(Color intensity, Ray ray, int depth, bool inside) {
		if (depth > D_MAX) return;

		Intersection firstIntersection = GetFirstIntersect(ray);
		int index = firstIntersection.GetIndex();
		if (index == -1) return;

		Object *object = objects[index];

		Vector N = firstIntersection.GetNormal();
		// IntersectPoint = x
		Vector point = firstIntersection.GetIntersectPoint();
		// RAy DIR = v
		Vector rayDir = ray.GetV().GetNormalized();
		if (N * rayDir > 0) N = N * -1.0;

		Color objectN;
		Color fresnel;
		Color kt;
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

	void ToneMap(float *image) {

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
				image[i] = image[i] * D / I;
				image[i + 1] = image[i + 1] * D / I;
				image[i + 2] = image[i + 2] * D / I;
			}
		}
	}

	void Render(float *image) {
		Ray ray;
		if (camera != NULL) {
			if (KAUSZTIKA) {
				if (FOTON_FAJLBOL) {

					cout << endl << endl << "Fotonterkepek fajlbol..." << endl << endl;

					for (int i = 0; i < numOfObjects; i++) {
						if (objects[i]->IsDiffuse()) {
							char fileName[200];
							sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", FOTONSZAM, i);
							//sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", fotonszam, i);
							objects[i]->ReadPhotonMapFromFile(fileName);
						}
					}

					cout << "Kesz." << endl << endl;

				}
				else {

					cout << "Fotonloves start..." << endl << endl;

					double fx, fy, fz;
					double length, intensity;
					Ray photonRay;
					Color lightColor;
					Vector lightPosition;
					Vector direction;

					Intersection testIntersection;

					LightSource *actual;

					for (int i = 0; i < numOfLights; i++) {
						cout << i + 1 << ". lampa fotonloves..." << endl;

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

						for (int j = 0; j<FOTONSZAM; j++) {
							//for(unsigned long j=0;j<fotonszam;j++) {
							//if ( !(j% (FOTONSZAM / 100 )) ) cout << (double)(j*100) / FOTONSZAM << "%..." << endl;

							if (!(j % (FOTONSZAM / 100))) {
								//cout << j << endl;
								percent = (unsigned long)(j * 100) / FOTONSZAM;
								cout << percent << "%";
								if (percent > 0) {
									actualtime = glutGet(GLUT_ELAPSED_TIME);		// program inditasa ota eltelt ido
									ellapsed = actualtime - starttime;
									remaining = (ellapsed / percent) * (100 - percent);
									cout << ", becsult hatralevo ido: " << remaining / 1000 << " s" << endl;
								}
								else cout << endl;
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

					cout << endl << endl << "Fotonloves kesz." << endl << endl;

					if (FOTON_FAJLBA) {

						cout << endl << endl << "Fotonterkepek fajlba..." << endl << endl;

						for (int i = 0; i < numOfObjects; i++) {
							if (objects[i]->IsDiffuse()) {
								char fileName[200];
								sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", FOTONSZAM, i);
								//sprintf_s(fileName, 50, "fotonmap-%d-%d.dat", fotonszam, i);
								objects[i]->WritePhotonMapToFile(fileName);
							}
						}
						cout << "Kesz." << endl << endl;
					}

				}
			}
			cout << "Sugarkovetes start..." << endl << endl;

			Color actual;
			for (int z = 0; z < KEPMAGASSAG; z++) {

				cout << (z + 1) << ". sor..." << endl;

				for (int x = 0; x < KEPSZELESSEG; x++) {
					bool willFail = false;
					ray = camera->GetRay(x, z);

					if (x == 260 && z == 175) {
						willFail = true;
					}

					actual = Trace(ray, 0, false);

					int actualIndex = (KEPMAGASSAG - 1 - z) * 3 * KEPSZELESSEG + x * 3;
					image[actualIndex] = actual.GetR();
					image[actualIndex + 1] = actual.GetG();
					image[actualIndex + 2] = actual.GetB();
				}
			}
		}
	}

};

// IMAGE
float finalImage[KEPSZELESSEG*KEPMAGASSAG * 3];

// COLORS
Color black(0.0, 0.0, 0.0);
Color darkGray(0.1f, 0.1f, 0.1f);
Color yellow(0.98f, 0.94f, 0.24f);
Color encyan(0.31f, 0.91f, 0.91f);
Color green(0.0, 0.1f, 0.0);
Color red(1.0, 0.0, 0.0);
Color white(1.0, 1.0, 1.0);
Color blue(0.0, 0.0, 1.0);
Color lightBlue(0.41f, 0.79f, 0.91f);
Color brown(0.78f, 0.58f, 0.47f);

// n and k values
Color copperN(0.2f, 1.1f, 1.2f);
Color copperK(3.6f, 2.6f, 2.3f);
Color silverN(0.14f, 0.16f, 0.13f);
Color silverK(4.1f, 2.3f, 3.1f);
Color goldN(0.17f, 0.35f, 1.5f);
Color goldK(3.1f, 2.7f, 1.9f);
Color transparentN(1.1f, 1.1f, 1.1f);
Color transparentK(0.0, 0.0, 0.0);
Color glassN(1.5f, 1.5f, 1.5f);
Color glassK(0.0, 0.0, 0.0);
Color diamondN(2.4f, 2.4f, 2.4f);
Color diamondK(0.0, 0.0, 0.0);

// LIGHT1
LightSource light1;
//Vector light1Position(800.0, 60.0, 500.0);
Vector light1Position(400.0, 150.0, 500.0);
double light1Intensity = 500000.0;
Color light1Color = white;

// OBJ1 = TABLE
Rect table;
PhotonMap tablePhotonMap;

Vector corner1(-600.0, 300.0);
Vector corner2(300.0, -600.0);
Vector corner3(1200.0, 300.0);
Vector corner4(300.0, 1200.0);

// OBJ2 = Golden ring
Cylinder goldenRing;

Vector goldenRingCenter(200.0, 100.0, 0.0);
double goldernRingRadius = 300.0;
double goldernRingHeight = 125.0;

// OBJ3 = Silver ring
Cylinder silverRing;

Vector silverRingCenter(450.0, 100.0, 0.0);
double silverRingRadius = 300.0;
double silverRingHeight = 125.0;

// OBJ4 = Copper ring
Cylinder copperRing;

Vector copperRingCenter(150.0, 100.0, 0.0);
double copperRingRadius = 300.0;
double copperRingHeight = 100.0;

// OBJ5 = Copper disk
Cylinder copperDisk;

Vector copperDiskCenter(0.0, 800.0, 0.0);
double copperDiskRadius = 200.0;
double copperDiskHeight = 75.0;

// OBJ6 = Diamond disk
Cylinder diamondDisk;

//Vector diamondDiskCenter(550.0, 100.0, 0.0);
Vector diamondDiskCenter(200.0, 150.0, 0.0);
double diamondDiskRadius = 150.0;
double diamondDiskHeight = 400.0;

// OBJ7 = Diamond1
Diamond diamond1;
Ellipsoid diamond1Bounding;
TriangleMesh diamond1TriangleMesh;
Triangle diamond1Triangles[MAX_HAROMSZOG_SZAM];

Vector diamond1Center(550.0, 0.0, 130.0);
double diamond1Width = 125.0;
double diamond1Depth = 125.0;
double diamond1Height = 100.0;

// OBJ8 = Diamond2
Diamond diamond2;
Ellipsoid diamond2Bounding;
TriangleMesh diamond2TriangleMesh;
Triangle diamond2Triangles[MAX_HAROMSZOG_SZAM];

Vector diamond2Center(0, 100.0, 115.0);
double diamond2Width = 125.0;
double diamond2Depth = 125.0;
double diamond2Height = 100.0;

// OBJ9 = Diamond Ellipsoid
Ellipsoid diamondElli;
Vector diamondElliCenter(200.0, 150.0, 200.0);
double diamondElliWidth = 75.0;
double diamondElliDepth = 75.0;
double diamondElliHeight = 200.0;

// OBJ10 = Diffuse Ellipsoid
Ellipsoid diffuseElli;
Vector diffuseElliCenter(200.0, 150.0, 200.0);
double diffuseElliWidth = 75.0;
double diffuseElliDepth = 75.0;
double diffuseElliHeight = 200.0;

// OBJ11 = Heart Shape
HeartShape heart;

// CAMERA
Camera camera;
//Vector lookAt(300.0,-300.0,300.0);
double r = 1.0;
double angle = 45;
double diffX = r * cos(angle * 180 / PI);
double diffY = r * sin(angle * 180 / PI);
//Vector lookAt(0.0 + diffX, 575.0 + diffY, 0.0);
Vector lookAt(0.0, 583.0, 0.0);
double planeWidth = 600.0;
double planeHeight = 600.0;
double planeAngle1 = 0.0;
//double planeAngle2 = 25.0;
double planeAngle2 = 0.0;
double fovDegree = 54.0;
//double fovDegree = 90.0;

// SCENE
Scene scene;
Color ambient = lightBlue;

void buildDiamondTriangleMesh(TriangleMesh *mesh, Triangle *triangles, Vector center, double width, double depth, double height) {
	// Diamond1

	mesh->Clear();

	int i = 0, u, v;
	double uStart, uEnd, vStart, vEnd;

	Vector vertex1, vertex2, vertex3, vertex4;

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
	//asztal.Diffuse().Kd(darkGray).Ka(darkGray * PI).Ks(white).Shine(1000.0).ProcMode(Object::ProceduralMode::MULTIPLY);
	table.Diffuse().Kd(darkGray).Ka(darkGray * PI).Ks(white).Shine(1000.0).ProcMode(Object::ProceduralMode::NONE);
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
	Vector right(halfPlaneWidth * cos1, halfPlaneWidth * sin1, 0.0);
	Vector up(0.0, halfPlaneHeight * sin2, halfPlaneHeight * cos2);

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

void writeImageToFile() {
	cout << "BMP file iras kezdete..." << endl << endl;

	ofstream file(FILE_NEV);

	if (file.is_open()) {

		int BMPheadersize = 14;
		int DIBheadersize = 40;

		// Egy pixel 3 bájt
		int RowSize = KEPSZELESSEG * 3;

		int RowPadding = 4 - ((KEPSZELESSEG * 3) % 4);

		if (RowPadding == 4) RowPadding = 0;

		// A teljes kép mérete bájtokban
		int BitmapSize = KEPMAGASSAG * (RowSize + RowPadding);

		int FileSize = BMPheadersize + DIBheadersize + BitmapSize;

		/// BMP HEADER KEZDETE

		// BMP formátumot jelölõ karaktersorozat
		file << (unsigned char)'B';
		file << (unsigned char)'M';

		// BMP fájl mérete little-endian formátumban:
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(FileSize % 256);
			FileSize /= 256;
		}

		// Reserved bitek helyére 0-k
		for (int i = 0; i < 4; i++) file << (unsigned char)0;

		int BitmapOffset = BMPheadersize + DIBheadersize;

		// Kép adatainak offszete little-endian formátumban
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(BitmapOffset % 256);
			BitmapOffset /= 256;
		}

		/// DIB HEADER KEZDETE

		// DIB header mérete little-endian-ban
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(DIBheadersize % 256);
			DIBheadersize /= 256;
		}

		int Width = KEPSZELESSEG;
		int Height = KEPMAGASSAG;

		// Képszélesség little-endian-ban
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(Width % 256);
			Width /= 256;
		}

		// Képmagasság little-endian-ban
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(Height % 256);
			Height /= 256;
		}

		// Color-plane-ek száma = 1  (  01 00 =>  00 01)
		file << (unsigned char)1;
		file << (unsigned char)0;

		// Pixelenkénti bitek száma
		file << (unsigned char)24;
		file << (unsigned char)0;

		// Sima RGB formátum
		for (int i = 0; i < 4; i++) file << (unsigned char)0;

		// Képadat mérete
		for (int i = 0; i < 4; i++) {
			file << (unsigned char)(BitmapSize % 256);
			BitmapSize /= 256;
		}

		// VízColortes felbontás
		file << (unsigned char)13;
		file << (unsigned char)11;
		file << (unsigned char)0;
		file << (unsigned char)0;

		// Függõleges felbontás
		file << (unsigned char)13;
		file << (unsigned char)11;
		file << (unsigned char)0;
		file << (unsigned char)0;

		// Színek szám a palettán
		for (int i = 0; i < 4; i++) file << (unsigned char)0;

		// Fontos színek szám a palettán
		for (int i = 0; i < 4; i++) file << (unsigned char)0;

		/// NYERS ADAT

		int index = 0;

		double tempR;
		double tempG;
		double tempB;

		int R, G, B;

		for (int i = 0; i < KEPMAGASSAG; i++) {
			for (int j = 0; j < KEPSZELESSEG; j++) {

				if (j == 540 && i == 361) {
					int a = 8;
				}
				index = i * 3 * KEPSZELESSEG + j * 3;

				tempR = finalImage[index] * 255;
				tempG = finalImage[index + 1] * 255;
				tempB = finalImage[index + 2] * 255;

				R = min((int)floor(tempR + 0.5f), 255);
				G = min((int)floor(tempG + 0.5f), 255);
				B = min((int)floor(tempB + 0.5f), 255);

				// B, G, R a sorrend
				file << (unsigned char)B;
				file << (unsigned char)G;
				file << (unsigned char)R;

			}

			// Ha van padding, akkor írunk nullákat

			for (int i = 0; i < RowPadding; i++) file << (unsigned char)0;
		}

		// Elvileg kész a BMP fájlunk
		file.close();
	}

	cout << "BMP file iras vege..." << endl << endl;
}

// Inicializacio, a program futasanak kezdeten, az OpenGL kontextus letrehozasa utan hivodik meg (ld. main() fv.)
void onInitialization() {

	initLights();
	initObjects();
	initCamera();
	initScene();

	scene.Render(finalImage);
	if (TONEMAP) scene.ToneMap(finalImage);

	// Fájlba írás
	if (FILE_IRAS) writeImageToFile();

	/*
	Polinom poli1(6);

	poli1.SetCoefficient(6,8);
	poli1.SetCoefficient(5,64);
	poli1.SetCoefficient(4,5);
	poli1.SetCoefficient(3,-56);
	poli1.SetCoefficient(2,132);
	poli1.SetCoefficient(1,-41);
	poli1.SetCoefficient(0,72);

	double x0 = -1;

	std::cout << "Polinom 1: " << poli1 << std::endl << std::endl;

	RealRoots roots = poli1.getRootsWithNewton(x0);

	std::cout << "Init guess = " << x0 << std::endl << std::endl;

	std::cout << "Roots: " << std::endl << roots << std::endl << std::endl;
	*/



	glutPostRedisplay();
}

void onDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(KEPSZELESSEG, KEPMAGASSAG, GL_RGB, GL_FLOAT, finalImage);
	glutSwapBuffers();
}

void onKeyboard(unsigned char key, int x, int y) {
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

// Eger esemenyeket lekezelo fuggveny
void onMouse(int button, int state, int x, int y) {
	//if (button == GLUT_LEFT && state == GLUT_DOWN);  // A GLUT_LEFT_BUTTON / GLUT_RIGHT_BUTTON illetve GLUT_DOWN / GLUT_UP
}

// `Idle' esemenykezelo, jelzi, hogy az ido telik, az Idle esemenyek frekvenciajara csak a 0 a garantalt minimalis ertek
void onIdle() {
	long time = glutGet(GLUT_ELAPSED_TIME);		// program inditasa ota eltelt ido

}

int main(int argc, char **argv) {

	Vector a;
	a.X(1.0).Y(2.0);

	/*
	if ( argc != 2 ) return 0;

	char *param = argv[1];

	int hossz;

	for ( hossz = 0; param[hossz] != 0; hossz++) if ( (param[hossz] > '9') || (param[hossz] < '0' ) ) return 0;

	int szorzo = 1;

	for ( hossz--, fotonszam = 0 ; hossz >= 0; hossz--, szorzo*=10 ) fotonszam += szorzo*(param[hossz]-'0');

	cout << fotonszam << "foton lesz... " << endl << endl;

	*/

	glutInit(&argc, argv); 				// GLUT inicializalasa
	glutInitWindowSize(KEPSZELESSEG, KEPMAGASSAG);			// Alkalmazas ablak kezdeti merete 600x600 pixel 
	glutInitWindowPosition(100, 100);			// Az elozo alkalmazas ablakhoz kepest hol tunik fel
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);	// 8 bites R,G,B,A + dupla buffer + melyseg buffer

	glutCreateWindow("Grafika hazi feladat");		// Alkalmazas ablak megszuletik es megjelenik a kepernyon

	glMatrixMode(GL_MODELVIEW);				// A MODELVIEW transzformaciot egysegmatrixra inicializaljuk
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);			// A PROJECTION transzformaciot egysegmatrixra inicializaljuk
	glLoadIdentity();

	onInitialization();					// Az altalad irt inicializalast lefuttatjuk

	glutDisplayFunc(onDisplay);				// Esemenykezelok regisztralasa
	glutMouseFunc(onMouse);
	glutIdleFunc(onIdle);
	glutKeyboardFunc(onKeyboard);

	glutMainLoop();					// Esemenykezelo hurok

	return 0;
}