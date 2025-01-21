#pragma once

#include "rectangularvectors.h"

namespace fea
{

	// Symbolic constants
	class FE
	{

		// Main method of application: JFEM/JMGEN/JVIS
	public:
		static int main;
		static constexpr int JFEM = 0;
			static constexpr int JMGEN = 1;
			static constexpr int JVIS = 2;

		static constexpr int maxNodesPerElem = 20;

		// Big value for displacements boundary conditions
		static double bigValue;
		// Solution tuning
		static bool tunedSolver;

		// Error tolerance for PCG solution of FE equation system
		static constexpr double epsPCG = 1.e-10;
		// Constants for PCG solution method
		static int maxRow2D;
			static int maxRow3D;
			static int maxIterPcg;

		// Integration scheme for elastic-plastic problems
		static bool epIntegrationTANGENT;
	};
}
