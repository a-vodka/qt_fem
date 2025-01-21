#pragma once

#include <vector>
#include "rectangularvectors.h"

namespace elem
{

	// Stresses and equivalent strains at integration point
	class StressContainer
	{

		// Accumulated stress
	public:
		std::vector<double> sStress;
		// StressContainer increment
		std::vector<double> dStress;
		// Accumulated equivalent plastic strain
		double sEpi = 0;
		// Equivalent plastic strain increment
		double dEpi = 0;

		StressContainer(int nDim);

	};
}
