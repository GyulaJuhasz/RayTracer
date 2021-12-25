#ifndef Color_H
#define Color_H

namespace raytracer {

	class Color {
		double r;
		double g;
		double b;

	public:

		Color(double r = 0.0, double g = 0.0, double b = 0.0);

		double GetR() const;
		double GetG() const;
		double GetB() const;

		Color& R(double newR);

		Color& G(double newG);

		Color& B(double newB);

		void Set(double cR = 0.0, double cG = 0.0, double cB = 0.0);

		Color& operator=(const Color& other);

		Color operator+(const Color& other);
		Color operator-(const Color& other);
		Color operator*(const double& c);
		Color operator/(const double& c);
		Color operator*(const Color& other);
		Color operator/(const Color& other);

	};

}

#endif // !Color_H