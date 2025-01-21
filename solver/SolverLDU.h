#pragma once

#include "Solver.h"
#include <vector>
#include "rectangularvectors.h"
#include "../fea/FE.h"

namespace solver
{

	using namespace fea;
	using namespace model;
	using namespace elem;


	// Profile LDU symmetric solver.
	// Upper symmetric part of the global stiffness matrix is
	// stored by columns of variable height (profile storage).
	class SolverLDU : public Solver
	{

		// Pointers to matrix columns
	private:
		std::vector<int> pcol;
		// Global stiffness matrix
		std::vector<double> A;

		// Constructor for LDU symmetric solver.
	public:
		SolverLDU();

		// Calculate profile of GSM: set column pointers pcol[]
	private:
		void setProfile();

		// Assemble element matrix to the global stiffness matrix
	public:
		void assembleESM() override;

		// Solve equation system by direct LDU method.
		// x - right-hand side/solution (in/out)
		int solve(std::vector<double> &x) override;

		// Apply displacement boundary conditions
	private:
		void displacementBC();

		// LDU decomposition for symmetric matrix
		int lduDecomposition();

		// Forward reduction and backsubstitution
		void lduFrwdBksb(std::vector<double> &x);

		// Tuned LDU decomposition for symmetric matrix
		// (block-block tuning, block size = 2 - 2D; 3 - 3D)
		int lduDecompositionTuned();

	};

}
