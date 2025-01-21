#pragma once

#include "Material.h"
#include <string>
#include <vector>
#include "rectangularvectors.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace elem { class Element; }

namespace material
{

	using Element = elem::Element;

	// Constitutive relations for elastic material
	class ElasticMaterial : public Material
	{

	public:
		static std::vector<double> deps;
		static std::vector<double> dsig;
		// Length of strain and stress vectors
		static int lv;

		ElasticMaterial(const std::wstring &stressState);

		// Hooke's law: increment of stress due to
		// increment of strain
		void strainToStress(Element *elm, int ip) override;

		// Compute elasticity matrix emat
		void elasticityMatrix(std::vector<std::vector<double>> &emat) override;

		// Elasticity 3D matrix emat [6][6]
		virtual void elasticityMatrix3D(std::vector<std::vector<double>> &emat);

		// Elasticity 2D matrix emat [4][4]
		virtual void elasticityMatrix2D(std::vector<std::vector<double>> &emat);

	};

}
