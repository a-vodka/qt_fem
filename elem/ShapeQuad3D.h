#pragma once

#include <vector>
#include "rectangularvectors.h"

namespace elem
{


	// Quadratic 3D shape functions and their derivatives
	class ShapeQuad3D
	{

		// Degeneration check.
		// The only degeneration is: 0=7=6,8=11,12=19=18
	public:
		static int degeneration(std::vector<int> &ind);

		// Shape functions.
		// xi, et, ze - local coordinates;
		// ind - element connectivities;
		// n - shape functions (out)
		static void shape(double xi, double et, double ze, std::vector<int> &ind, std::vector<double> &n);

		// Derivatives of shape functions
		//    with respect to global coordinates xy.
		// xi, et, ze - local coordinates;
		// ind - element connectivities;
		// xy - nodal coordinates;
		// dnxy - derivatives of shape functions (out);
		// returns  determinant of the Jacobian matrrix
		static double deriv(double xi, double et, double ze, std::vector<int> &ind, std::vector<std::vector<double>> &xy, std::vector<std::vector<double>> &dnxy);

		// Two-dimensional shape functions and derivatives
		//    for a face of 3d 8-20n element.
		// xi, et - local coordinates;
		// ind - element connectivities;
		// an - shape functions (out);
		// dn  derivatives of shape functions
		//      with respect to xi and et (out)
		static void shapeDerivFace(double xi, double et, std::vector<int> &ind, std::vector<double> &an, std::vector<std::vector<double>> &dn);

	};
}
