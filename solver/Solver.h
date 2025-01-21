#pragma once

#include <vector>
#include "rectangularvectors.h"
#include "../elem/Element.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace model { class FeModel; }

namespace solver
{

	using namespace elem;
	using namespace model;

	// Solution of the global equation system
	class Solver
	{

	public:
		static FeModel *fem;
		// number of equations
		static int neq;
		// length of global stiffness matrix
		static int lengthOfGSM;
		// elem connectivities - degrees of freedom
		std::vector<int> indf;
		// number of degrees of freedom for element
		int nindf = 0;
		// Indicator of new global matrix
		bool newMatrix = false;

	public:
		enum class Solvers
		{
//JAVA TO C++ CONVERTER TODO TASK: Enum value-specific class bodies are not converted by Java to C++ Converter:
//			ldu {Solver create() {return new SolverLDU();} },
//JAVA TO C++ CONVERTER TODO TASK: Enum value-specific class bodies are not converted by Java to C++ Converter:
//			pcg {Solver create() {return new SolverPCG();} };
//JAVA TO C++ CONVERTER TODO TASK: Enum methods are not converted by Java to C++ Converter:
//			abstract Solver create();
		};

	public:
		static Solvers solver;

		static Solver *newSolver(FeModel *fem);

		// Assemble global stiffnes matrix
		virtual void assembleGSM();

		// Add element stiffness matrix to GSM
		virtual void assembleESM();

		// Solve global equation system
		// x - right-hand side/solution (in/out)
		virtual int solve(std::vector<double> &x);

	};
}
