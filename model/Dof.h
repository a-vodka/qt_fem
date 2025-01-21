#pragma once

#include "rectangularvectors.h"

namespace model
{

	// Degree of freedom.
	class Dof
	{

	public:
		int dofNum = 0;
		double value = 0;

		Dof(int dofNum, double value);

	};

}
