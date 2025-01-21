#include "FeStress.h"
#include "FeModel.h"
#include "FeLoad.h"
#include "../elem/Element.h"
#include "../material/Material.h"
#include "Dof.h"
#include "elem/StressContainer.h"

//#include "../util/FePrintWriter.h"


namespace model
{
	using namespace elem;
	using namespace material;
	using namespace util;
double FeStress::relResidNorm = 0;

	FeStress::FeStress(FeModel *fem)
	{
		this->fem = fem;
	}

	void FeStress::computeIncrement()
	{

		// Accumulate solution vector in displacement increment
		for (int i = 0; i < fem->nEq; i++)
		{
			FeLoad::dDispl[i] += FeLoad::RHS[i];
		}

		// Compute stresses at reduced integration points
		for (int iel = 0; iel < fem->nEl; iel++)
		{
			Element *elm = fem->elems[iel];
			elm->setElemXyT();
			elm->disAssembleElemVector(FeLoad::dDispl);

            for (size_t ip = 0; ip < elm->str.size(); ip++)
			{
                Material* mat = fem->materials[elm->matName];
                mat->strainToStress(elm, ip);
			}
		}
	}

	bool FeStress::equilibrium(int iter)
	{

		if (fem->physLaw == FeModel::PhysLaws::elastic || iter == FeLoad::maxIterNumber)
		{
			return true;
		}
		// Assemble residual vector to right-hand side
		for (int i = 0; i < fem->nEq; i++)
		{
			FeLoad::RHS[i] = FeLoad::spLoad[i] + FeLoad::dpLoad[i];
		}
		Element *elm;
		for (int iel = 0; iel < fem->nEl; iel++)
		{
			elm = fem->elems[iel];
			elm->setElemXy();
			elm->equivStressVector();
			elm->assembleElemVector(Element::evec,FeLoad::RHS);
		}
        // Displacement boundary conditions
        auto it = fem->defDs.begin();
        while (it != fem->defDs.end())
        {
            auto d = it;
			FeLoad::RHS[d->dofNum - 1] = 0;
			it++;
		}
		// Relative residual norm
		double dpLoadNorm = vectorNorm(FeLoad::dpLoad);
		if (dpLoadNorm < 1e-30)
		{
			dpLoadNorm = vectorNorm(FeLoad::dhLoad);
		}
		relResidNorm = vectorNorm(FeLoad::RHS) / dpLoadNorm;
		return relResidNorm < FeLoad::residTolerance;
	}

	double FeStress::vectorNorm(std::vector<double> &v)
	{

		double norm = 0;
		for (auto aV : v)
		{
			norm += aV * aV;
		}
		return std::sqrt(norm);
	}

	void FeStress::accumulate()
	{

		for (int i = 0; i < fem->nEq; i++)
		{
			FeLoad::spLoad[i] += FeLoad::dpLoad[i];
			FeLoad::sDispl[i] += FeLoad::dDispl[i];
		}
		for (int iel = 0; iel < fem->nEl; iel++)
		{
			 fem->elems[iel]->accumulateStress();
		}
	}

	void FeStress::writeResults()
	{

        if (fem->nDim == 2)
        {
             std::wcout << L" Node             ux             uy\n";
        }
        else
        {
             std::wcout << L" Node             ux             uy             uz\n";
        }
        for (int i = 0; i < fem->nNod; i++)
        {
             std::wcout << i + 1<<"\t";
             for (int j = 0; j < fem->nDim; j++)
             {
                std::wcout << FeLoad::sDispl[fem->nDim * i + j] << "\t";
             }
             std::wcout<<std::endl;
        }

        std::wcout <<L"\n\nStresses\n";
        for (int iel = 0; iel < fem->nEl; iel++)
        {
             if (fem->nDim == 2)
             {
                std::wcout << L"\nEl " << iel + 1 <<"\t sxx\t syy\t sxy\t szz\t epi";
             }
             else
             {
                std::wcout << L"\nEl " << iel + 1 << L"\t sxx\t syy\t szz\t sxy\t syz\t szx\t epi";
             }
             for (auto aStr : fem->elems[iel]->str)
             {
                std::wcout <<L"\n\t";
                for (int i = 0; i < 2 * fem->nDim; i++)
                {
                    std::wcout << aStr->sStress[i] <<'\t';
                }
                std::wcout << aStr->sEpi;
             }
        }
        std::wcout <<L"\n";

/*
		std::wstring fileResult = fea::Jfem::fileOut + L"." + FeLoad::loadStepName;
		FePrintWriter tempVar();
		PrintWriter *PR = (&tempVar)->getPrinter(fileResult);

		PR->printf(L"Displacements\n\n");
		if (fem->nDim == 2)
		{
			PR->printf(L" Node             ux             uy");
		}
		else
		{
			PR->printf(L" Node             ux             uy" + L"             uz");
		}
		for (int i = 0; i < fem->nNod; i++)
		{
			PR->printf(L"\n%5d", i + 1);
			for (int j = 0; j < fem->nDim; j++)
			{
			  PR->printf(L"%15.6e", FeLoad::sDispl[fem->nDim * i + j]);
			}
		}

		PR->printf(L"\n\nStresses\n");
		for (int iel = 0; iel < fem->nEl; iel++)
		{
			if (fem->nDim == 2)
			{
				PR->printf(L"\nEl %4d     sxx            syy" + L"            sxy            szz" + L"            epi", iel + 1);
			}
			else
			{
				PR->printf(L"\nEl %4d     sxx            syy" + L"            szz            sxy" + L"            syz            szx" + L"            epi", iel + 1);
			}
			for (auto aStr : *fem->elems[iel]->str)
			{
				PR->printf(L"\n");
				for (int i = 0; i < 2 * fem->nDim; i++)
				{
					PR->printf(L"%15.8f", aStr->sStress[i]);
				}
				PR->printf(L"%15.8f", aStr->sEpi);
			}
		}
		PR->close();
*/
	}

	void FeStress::readResults(const std::wstring &resultFile, std::vector<double> &displ)
	{
/*
		if (resultFile == L"")
		{
			return;
		}

		FeScanner *RD = new FeScanner(resultFile);
		// Read displacements
		RD->moveAfterLineWithWord(L"node");
		for (int i = 0; i < fem->nNod; i++)
		{
			RD->readInt();
			for (int j = 0; j < fem->nDim; j++)
			{
				displ[fem->nDim * i + j] = RD->readDouble();
			}
		}
		// Read stresses
		for (int iel = 0; iel < fem->nEl; iel++)
		{
			RD->moveAfterLineWithWord(L"el");
			for (auto aStr : *fem->elems[iel]->str)
			{
				for (int i = 0; i < 2 * fem->nDim; i++)
				{
					aStr->sStress[i] = RD->readDouble();
				}
				aStr->sEpi = RD->readDouble();
			}
		}
		RD->close();

		delete RD;
*/
	}

}
