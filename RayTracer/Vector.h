#ifndef VECTOR_H
#define VECTOR_H

namespace raytracer {

	class Vector {
		double x;
		double y;
		double z;
		double h;

	public:
		Vector(double x = 0.0, double y = 0.0, double z = 0.0, double h = 1.0);

		double GetX() const;
		double GetY() const;
		double GetZ() const;
		double GetH() const;

		double GetRawX() const;
		double GetRawY() const;
		double GetRawZ() const;

		Vector& X(double newX);
		Vector& Y(double newY);
		Vector& Z(double newZ);
		Vector& H(double newH);

		void Set(double vX = 0.0, double vY = 0.0, double vZ = 0.0, double vH = 1.0);

		double GetLenght();

		Vector operator+(const Vector& v);
		Vector operator-(const Vector& v);
		Vector operator%(const Vector& v);
		Vector operator*(const double& c);
		double operator*(const Vector& v);
		Vector operator/(const double& c);

		Vector& operator=(const Vector& v);

		void Normalize();

		Vector GetNormalized();
	};

}

#endif // !VECTOR_H