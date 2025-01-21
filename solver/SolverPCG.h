#pragma once

#include "Solver.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include "rectangularvectors.h"
#include "../fea/FE.h"

namespace solver
{

	using namespace fea;
	using namespace model;
	using namespace elem;


	// Preconditioned conjugate gradient (PCG) solver.
	// Matrix storage: sparse-row format
	class SolverPCG : public Solver
	{

		// Pointers to matrix rows
	private:
		std::vector<int> prow;
		// Column numbers for non-zero values in A
		std::vector<int> coln;
		// Nonzero values of the global stiffness matrix by rows
		std::vector<double> A;
		// Working arrays
		std::vector<double> b; std::vector<double> r; std::vector<double> w; std::vector<double> p; std::vector<double> md;

		// Constructor for PCG solver.
	public:
		SolverPCG();

		// Set sparse row structure for storage of nonzero
		//  coefficients of the global stiffness matrix.
	private:
		void setSparseRowStructure();

		// Assemble element matrix to the global stiffness matrix
	public:
		void assembleESM() override;

		// Solve equation system by PCG method.
		// x - right-hand side/solution (in/out)
		int solve(std::vector<double> &x) override;

		// Apply displacement boundary conditions:
	private:
		void displacementBC();

		// PCG solution method.
		// x - right-hand side/solution (in/out)
	public:
		virtual int pcg(std::vector<double> &x);

		// Diagonal preconditioner md = (Diag(A))-1.
	private:
		void diagonalPreconditioner();

		//  Sparse matrix-vector product y = A*x.
		void matrixVectorProduct(std::vector<double> &x, std::vector<double> &y);

	};

}
