#include "../solver/Solver.h"
#include "../model/FeModel.h"
#include "../elem/Element.h"
#include "../fea/FE.h"
#include "solver/SolverLDU.h"

namespace solver
{
	using namespace elem;
	using namespace model;
	using FE = fea::FE;
    FeModel *Solver::fem;
    int Solver::neq = 0;
    int Solver::lengthOfGSM = 0;
    //Solvers Solver::solver = Solvers::ldu;

	Solver *Solver::newSolver(FeModel *fem)
	{
		Solver::fem = fem;
        neq = fem->nEq;
        return new SolverLDU();
	}

	void Solver::assembleGSM()
	{
		Element *elm;
		indf = std::vector<int>(FE::maxNodesPerElem * fem->nDf);

		for (int iel = 0; iel < fem->nEl; iel++)
		{
            for (size_t i = 0; i < fem->elems[iel]->ind.size(); i++)
			{
				for (int k = 0; k < fem->nDf; k++)
				{
					indf[fem->nDf * i + k] = (fem->elems[iel]->ind[i] - 1) * fem->nDf + k + 1;
				}
			}
            nindf = fem->elems[iel]->ind.size() * fem->nDf;
			elm = fem->elems[iel];
			elm->setElemXy();
			elm->stiffnessMatrix();
			assembleESM();
		}
		// Indicate that new global matrix appeared
		newMatrix = true;
	}

	void Solver::assembleESM()
	{
	}

	int Solver::solve(std::vector<double> &x)
	{
		return 0;
	}
}
