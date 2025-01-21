#include "StressContainer.h"

namespace elem
{

	StressContainer::StressContainer(int nDim)
	{
		sStress = std::vector<double>(2 * nDim);
		dStress = std::vector<double>(2 * nDim);
	}
}
