#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "rectangularvectors.h"
#include "../elem/Element.h"
#include "../util/UTIL.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace model { class FeModel; }

namespace model
{

	using namespace elem;
	using namespace material;
	using namespace util;


	// Stress increment due to displacement increment
	class FeStress
	{

	public:
		static double relResidNorm;
	private:
		FeModel *fem;

		// Constructor for stress increment.
		// fem - finite element model
	public:
		virtual ~FeStress()
		{
			delete fem;
		}

		FeStress(FeModel *fem);

		// Compute stress increment for the finite element model
		virtual void computeIncrement();

		// Check equilibrium and assemble residual vector.
		// iter - number of iterations performed
		virtual bool equilibrium(int iter);

		// Returns norm of a vector v
		virtual double vectorNorm(std::vector<double> &v);

		// Accumulate loads, temperature and stresses
		virtual void accumulate();

		// Write results to a file.
		virtual void writeResults();

		// Read results from a file.
		// displ - displacements for the finite element model (out)
		virtual void readResults(const std::wstring &resultFile, std::vector<double> &displ);

	};

}
