#pragma once

#include "ElasticMaterial.h"
#include <string>
#include <vector>
#include <cmath>
#include "rectangularvectors.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:

namespace elem { class Element; }

namespace material
{

    using Element = elem::Element;
    class Midpoint;

	//  Constitutive relations for elastic-plastic material
	class ElasticPlasticMaterial : public ElasticMaterial
	{

		// Stress at the beginning of increment
	private:
		static std::vector<double> sig0;
		// Stress at the end of increment
		static std::vector<double> sig;
		// Derivatives of yield function
		static std::vector<double> a;
		// Elasticity matrix
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double[][] Emat = new double[6][6];
		static std::vector<std::vector<double>> Emat;
		// Ea = Emat*a
		static std::vector<double> Ea;
		// Plastic strain increment
		static std::vector<double> depsp;
		// Equivalent plastic strain
		static double epi;
		// Shear modulus
		static double G;
		static constexpr double beta = 0.1;
		Midpoint *midpoint;

	public:
		virtual ~ElasticPlasticMaterial()
		{
			delete midpoint;
		}

		ElasticPlasticMaterial(const std::wstring &stressState);

		// Elastic-plastic stress increment.
		// elm - element,
		// ip - integration point within element
		void strainToStress(Element *elm, int ip) override;

		// Compute elastic-plastic increment by tangent method.
		// Update stresses sig and equivalent plastic strain epi
	private:
		void tangentStressIncrement();

		// Yield function.
		// s - stresses,
		// epi - equivalent plastic strain,
		// returns  yield function value
		double yieldFunction(std::vector<double> &s, double epi);

		// Radius of yield surface Y = sY + k*ep^m
		double yieldRadius(double ep);

		// Derivatives of yield function.
		// s - stresses (in),
		// a - derivatives of yield function (out)
		void derivYieldFunc(std::vector<double> &s, std::vector<double> &a);

		// Returns slope of deformation curve
		// epi - equivalent plastic strain
		double slopeH(double epi);

		// Returns equivalent plastic strain
		// dp - pastic strains (in)
		double eqPlastStrain(std::vector<double> &dp);

		// Midpoint method for integration of constitutive
		// relations for the von Mises hardening material
	public:
		class Midpoint
		{
		private:
			ElasticPlasticMaterial *outerInstance;

		public:
			virtual ~Midpoint()
			{
				delete outerInstance;
			}

			Midpoint(ElasticPlasticMaterial *outerInstance);

			// Strain deviator
            double ed[6];
			// Stress deviator
            double sd[6];
			// Trial stress deviator
            double sdtr[6];
			// Yield function derivatives
            double a0[6];

            double sal[6];
            double salbar[6];
            double b[6];
			double alpha = 0.5;

			// Elastic-plastic stress increment by midpoint method.
			// Update stresses sig and
			//    equivalent plastic strain epi
			virtual void stressIncrement();

			// Compute deviator.
			// s - stress,
			// d - deviator (out),
			// returns  mean value
            virtual double deviator(double s[], double d[]);
            double deviator(std::vector<double> &s, double d[]);

			// Compute stress s = d + sm.
			// d - deviator,
			// sm - mean stress,
			// s - stress (out)
            virtual void stressFromDeviator(double d[], double sm, double s[]);

			// Returns norm = sqrt(Sij*Sij)
            virtual double norm(double s[]);

			// Returns dyadic product = aij*bij
            virtual double dyadicProduct(double a[], double b[]);

		};

	};

}
