#ifndef POINT_H
#define POINT_H

#include <math.h>

namespace estimator {

	class Point {
	public:
		Point(float x, float y, float z) : x(x), y(y), z(z) {}

		float x, y, z;

		void flipSign() {
			x *= -1;
			y *= -1;
			z *= -1;
		}

		void add(estimator::Point p) {
			x += p.x;
			y += p.y;
			z += p.z;
		}

		void multiply(int scalar) {
			x *= scalar;
			y *= scalar;
			z *= scalar;
		}

		float dot(estimator::Point p) {
			return (x * p.x) + (y * p.y) + (z * p.z);
		}

		void normalize() {
			float len = sqrt((x * x) + (y * y) + (z * z));
			x /= len;
			y /= len;
			z /= len;
		}
	};

}

#endif