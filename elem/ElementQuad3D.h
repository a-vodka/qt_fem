#pragma once

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

	// 3D 8-20 node isoparametric brick-type element.
	class ElementQuad3D : public Element
	{

	private:
		static std::vector<std::vector<int>> faceInd;
		static std::vector<double> an;
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double[][] dnxy = new double[20][3];
		static std::vector<std::vector<double>> dnxy;
		// Gauss rules for stiffness matrix, thermal vector,
		// surface load and stress integration
		static GaussRule *gk;
		static GaussRule *gh;
		static GaussRule *gf;
		static GaussRule *gs;

		// Constructor for 3D 20 node element.
	public:
		ElementQuad3D();

		// Compute stiffness matrix
		void stiffnessMatrix() override;

		// Compute thermal vector
		void thermalVector() override;

		// Set nodal equivalent of distributed face load to evec.
		// surLd - object describing element face load;
		// returns loaded element face
		// or -1 (loaded nodes does not match element face)
		int equivFaceLoad(ElemFaceLoad *surLd) override;

		// Compute equivalent stress vector (with negative sign)
		void equivStressVector() override;

		// Extrapolate stresses from integration points to nodes.
		// fip [8][6] - stresses at integration points;
		// fn [20][6] - stresses at nodes (out)
		void extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn) override;

		// Get local node numbers for element faces.
		// returns elementFaces[nFaces][nNodesOnFace]
		std::vector<std::vector<int>> getElemFaces() override;

		// Get strains at integration point.
		// ip - integration point number (stress);
		// returns  strain vector (ex, ey, ez, gxy, gyz, gzx)
		std::vector<double> getStrainsAtIntPoint(int ip) override;

		// Returns temperature at integration point (stress)
		double getTemperatureAtIntPoint(int ip) override;

	};
}
