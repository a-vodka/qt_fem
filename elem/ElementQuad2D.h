#pragma once

#define _USE_MATH_DEFINES
#include "Element.h"
#include <vector>
#include <cmath>
#include "rectangularvectors.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace util { class GaussRule; }
namespace model { class ElemFaceLoad; }

namespace elem
{

	using namespace model;
	using namespace material;
	using namespace util;

	// 2D quadratic isoparametric element (4-8 nodes)
	class ElementQuad2D : public Element
	{
		// Element edges (local numbers)
	private:
		static std::vector<std::vector<int>> faceInd;
		// Shape functions
		static std::vector<double> an;
		// Derivatives of shape functions
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double[][] dnxy = new double[8][2];
		static std::vector<std::vector<double>> dnxy;
		// Displacements differentiation matrix
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double[][] bmat = new double[4][16];
		static std::vector<std::vector<double>> bmat;
		// Elasticity matrix
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double[][] emat = new double[4][4];
		static std::vector<std::vector<double>> emat;
		// Thermal strains
		static std::vector<double> ept;
		// Radius in the axisymmetric problem
		static double r;
		// Gauss rules for stiffness matrix, thermal vector,
		// surface load and stress integration
		static GaussRule *gk;
		static GaussRule *gh;
		static GaussRule *gf;
		static GaussRule *gs;

		// Constructor for 2D quadratic element
	public:
		ElementQuad2D();

		// Compute stiffness matrix
		void stiffnessMatrix() override;

		// Set displacement differentiation matrix bmat.
		// xi, et - local coordinates,
		// returns  determinant of Jacobian matrix
	private:
		double setBmatrix(double xi, double et);

		// Compute thermal vector
	public:
		void thermalVector() override;

		// Set nodal equivalent of distributed face load to evec.
		// surLd - object describing element face load;
		// returns loaded element face
		// or -1 (loaded nodes do not match elem face)
		int equivFaceLoad(ElemFaceLoad *surLd) override;

		// Compute equivalent stress vector (with negative sign)
		void equivStressVector() override;

		// Extrapolate values from integration points to nodes.
		// fip [4][4] - values at integration points;
		// fn [8][4] - values at nodes (out)
		void extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn) override;

		// Get local node numbers for element faces.
		// returns elementFaces[nFaces][nNodesOnFace]
		std::vector<std::vector<int>> getElemFaces() override;

		// Get strains at integration point.
		// ip - integration point number (stress);
		// returns  strain vector (ex, ey, gxy, ez)
		std::vector<double> getStrainsAtIntPoint(int ip) override;

		// Get temperature at integration point (stress)
		double getTemperatureAtIntPoint(int ip) override;

	};

}
