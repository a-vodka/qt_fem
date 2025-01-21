#pragma once

#include <string>
#include <vector>
#include "rectangularvectors.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace elem { class Element; }

namespace material
{

	using Element = elem::Element;

	// Material constitutive relations
	class Material
	{

		// StressContainer state (plstrain/plstress/axisym/threed)
	public:
		std::wstring stressState;
		// Elasticity modulus
		double e = 0;
		// Poisson's ratio
		double nu = 0;
		// Thermal expansion
		double alpha = 0;
		// Yield stress
		double sY = 0;
		// Hardening coefficient
		double km = 0;
		// Hardening power
		double mm = 0;

        static Material *newMaterial(const wchar_t *matPhysLaw, const wchar_t *stressState);

		// Given strain increment at integration point ip
		// element elm, compute stress dsig increment
		virtual void strainToStress(Element *elm, int ip);

		// Set elastic properties
		virtual void setElasticProp(double e, double nu, double alpha);

		// Set plastic properties
		virtual void setPlasticProp(double sY, double km, double mm);

		// Returns Lame constant lambda
		virtual double getLambda();

		// Returns shear modulus
		virtual double getMu();

		// Returns Poisson's ratio
		virtual double getNu();

		// Returns thermal expansion coefficient
		virtual double getAlpha();

		// Compute elasticity matrix emat
		virtual void elasticityMatrix(std::vector<std::vector<double>> &emat);

	};

}
