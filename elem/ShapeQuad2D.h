#pragma once

#include <vector>
#include "rectangularvectors.h"

namespace elem
{


	// Quadratic 2D shape functions and their derivatives
	class ShapeQuad2D
	{

		// Degeneration check.
		// If element is triangular then the method returns
		// a local number (starting from 0) of the midside node
		// opposite to degenerated side.
		// ind - connectivity numbers
	public:
		static int degeneration(std::vector<int> &ind);

		// Shape functions.
		// xi, et - local coordinates;
		// ind[8] - element connectivities;
		// an[8] - shape functions (out)
		static void shape(double xi, double et, std::vector<int> &ind, std::vector<double> &an);

		// Derivatives of shape functions
		// with respect to global coordinates x and y.
		// xi, et - local coordinates;
		// ind[8] - element connectivities;
		// xy[8][2] - nodal coordinates;
		// dnxy[8][2] - derivatives of shape functions (out);
		// returns  determinant of the Jacobian matrrix
		 static double deriv(double xi, double et, std::vector<int> &ind, std::vector<std::vector<double>> &xy, std::vector<std::vector<double>> &dnxy);

		// One-dimensional quadratic shape functions and
		//    their derivatives in local coordinates
		// xi - local coordinate;
		// kmid - index of midside node (=0 no midside node);
		// an[3] - shape functions (out);
		// dndxi[3] - derivatives of shape functions (out)
		static void shapeDerivFace(double xi, int kmid, std::vector<double> &an, std::vector<double> &dndxi);

	};
}
