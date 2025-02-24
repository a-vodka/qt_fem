#include "SolverPCG.h"
#include "../fea/FE.h"
#include "../util/UTIL.h"
#include "../elem/Element.h"
#include "../model/Dof.h"
#include "../model/FeModel.h"
namespace solver
{
	using namespace fea;
	using namespace model;
	using namespace elem;
	using UTIL = util::UTIL;

	SolverPCG::SolverPCG()
	{

		// Set structure of matrix A
		setSparseRowStructure();

		A = std::vector<double>(prow[neq]);
		b = std::vector<double>(neq);
		r = std::vector<double>(neq);
		w = std::vector<double>(neq);
		p = std::vector<double>(neq);
		md = std::vector<double>(neq);
	}

	void SolverPCG::setSparseRowStructure()
	{
		prow = std::vector<int>(neq + 1);
		int lrow = fem->nDf * ((fem->nDf == 2) ? FE::maxRow2D : FE::maxRow3D);
		coln = std::vector<int>(neq * lrow);
		for (int i = 0; i <= fem->nNod; i++)
		{
			prow[i] = i * lrow * fem->nDf;
		}

		// Create nodal sparse-row matrix structure
		for (int i = 0; i < prow[fem->nNod]; i++)
		{
			coln[i] = -1;
		}
		// Diagonal entry - first in row
		for (int i = 0; i < fem->nNod; i++)
		{
			coln[prow[i]] = i;
		}
		for (int iel = 0; iel < fem->nEl; iel++)
		{
            for (auto anInd : fem->elems[iel]->ind)
			{
				if (anInd == 0)
				{
					continue;
				}
				int ii = anInd - 1; // Hyperrow
                for (auto anInd1 : fem->elems[iel]->ind)
				{
					if (anInd1 == 0)
					{
						continue;
					}
					int jj = anInd1 - 1; // Hypercolumn
					int k;
					for (k = prow[ii]; k < prow[ii + 1]; k++)
					{
						// If column already exists
						if (coln[k] == jj)
						{
							break;
						}
						if (coln[k] == -1)
						{
							coln[k] = jj;
							break;
						}
					}
					if (k == prow[ii + 1])
					{
                        UTIL::errorMsg(L"PCG sparse-row structure: not enough space for node " + std::to_wstring(ii));
					}
				}
			}
		}

		// Compress
		int p = 0;
		for (int i = 0; i < fem->nNod; i++)
		{
			int k = prow[i];
			prow[i] = p;
			for (int j = k; j < prow[i + 1]; j++)
			{
				if (coln[j] == -1)
				{
					break;
				}
				coln[p++] = coln[j];
			}
		}
		prow[fem->nNod] = p;

		// Transform to degrees of freedom
		int pdof = p * fem->nDf * fem->nDf;
		for (int i = fem->nNod - 1; i >= 0; i--)
		{
			int deln = (prow[i + 1] - prow[i]);
			p -= deln;
			for (int k = fem->nDf; k > 0; k--)
			{
				prow[i * fem->nDf + k] = pdof;
				pdof -= deln * fem->nDf;
				for (int j = prow[i + 1] - prow[i] - 1; j >= 0; j--)
				{
					for (int m = fem->nDf - 1; m >= 0; m--)
					{
						coln[pdof + j * fem->nDf + m] = coln[p + j] * fem->nDf + m;
					}
				}
			}
		}
		lengthOfGSM = static_cast<int>(prow[neq] * 1.5);
	}

	void SolverPCG::assembleESM()
	{
		for (int j = 0; j < nindf; j++)
		{
			int jj = indf[j] - 1;
			if (jj >= 0)
			{
				for (int i = 0; i < nindf; i++)
				{
					int ii = indf[i] - 1;
					if (ii >= 0)
					{
						// Sparse row format (full matrix)
						int k;
						for (k = prow[jj]; k < prow[jj + 1]; k++)
						{
							if (coln[k] == ii)
							{
								break;
							}
						}
						if (i <= j)
						{
							A[k] += Element::kmat[i][j];
						}
						else
						{
							A[k] += Element::kmat[j][i];
						}
					}
				}
			}
		}
	}

	int SolverPCG::solve(std::vector<double> &x)
	{

		if (newMatrix)
		{
			displacementBC();
			newMatrix = false;
		}
		return pcg(x);
	}

	void SolverPCG::displacementBC()
	{
        auto it = fem->defDs.begin();

        while (it!=fem->defDs.end())
		{
            auto d = it;
			int j = d->dofNum - 1;
			for (int k = prow[j]; k < prow[j + 1]; k++)
			{
				if (coln[k] == j)
				{
					A[k] = FE::bigValue;
				}
				else
				{
					A[k] = 0.0;
				}
			}
			it++;
		}
	}

	int SolverPCG::pcg(std::vector<double> &x)
	{

		diagonalPreconditioner();

		// Save x[] in b[] and calculate initial x
		for (int i = 0; i < neq; i++)
		{
			b[i] = x[i];
			x[i] = x[i] * md[i];
		}

		// r = b - A*x and initinal error
		matrixVectorProduct(x, r);
		for (int i = 0; i < neq; i++)
		{
			r[i] = b[i] - r[i];
		}
		double gamma0 = 1;
		double gammai = 1;
		int iter;
		for (iter = 0; iter < FE::maxIterPcg; iter++)
		{
			//  w = (M-1)*r
			for (int i = 0; i < neq; i++)
			{
				w[i] = md[i] * r[i];
			}
			// gam = (r,w)
			double gammai1 = gammai;
			gammai = 0;
			for (int i = 0; i < neq; i++)
			{
				gammai += r[i] * w[i];
			}
			if (iter == 0)
			{
				gamma0 = gammai;
				std::copy_n(w.begin(), neq, p.begin());
			}
			else
			{
				double rg = gammai / gammai1;
				for (int i = 0; i < neq; i++)
				{
					p[i] = w[i] + rg * p[i];
				}
			}
			// w = A*p
			matrixVectorProduct(p, w);
			double beta = 0;
			for (int i = 0; i < neq; i++)
			{
				beta += p[i] * w[i];
			}
			double alpha = gammai / beta;
			// Update x and r, calculate error
			for (int i = 0; i < neq; i++)
			{
				x[i] += alpha * p[i];
				r[i] -= alpha * w[i];
			}
			double err = std::sqrt(gammai / gamma0);
			if (err < FE::epsPCG)
			{
				return (iter + 1);
			}
		}
		return (iter);
	}

	void SolverPCG::diagonalPreconditioner()
	{
		for (int j = 0; j < neq; j++)
		{
			int i;
			for (i = prow[j]; i < prow[j + 1]; i++)
			{
					if (coln[i] == j)
					{
						break;
					}
			}
			md[j] = 1.0 / A[i];
		}
	}

	void SolverPCG::matrixVectorProduct(std::vector<double> &x, std::vector<double> &y)
	{
		if (FE::tunedSolver)
		{
			if (fem->nDf == 2)
			{ // tuned for nDf = 2
				for (int j = 0; j < neq; j++)
				{
					double s = 0;
					for (int i = prow[j]; i < prow[j + 1]; i += 2)
					{
						s += A[i] * x[coln[i]] + A[i + 1] * x[coln[i + 1]];
					}
					y[j] = s;
				}
			}
			else
			{ // tuned for nDf = 3
				for (int j = 0; j < neq; j++)
				{
					double s = 0;
					for (int i = prow[j]; i < prow[j + 1]; i += 3)
					{
						s += A[i] * x[coln[i]] + A[i + 1] * x[coln[i + 1]] + A[i + 2] * x[coln[i + 2]];
					}
					y[j] = s;
				}
			}
		}
		else
		{ // not tuned
			for (int j = 0; j < neq; j++)
			{
				double s = 0;
				for (int i = prow[j]; i < prow[j + 1]; i++)
				{
					s += A[i] * x[coln[i]];
				}
				y[j] = s;
			}
		}
	}
}
