#include "FeLoad.h"
//#include "FeModel.h"
#include "../util/FeScanner.h"
#include "../util/UTIL.h"
#include "ElemFaceLoad.h"
#include "../elem/Element.h"
#include "Dof.h"
#include "../fea/FE.h"

namespace model
{
	using namespace fea;
	using namespace elem;
	using namespace util;
FeModel *FeLoad::fem;

	FeLoad::FeLoad(FeModel *fem)
	{
		FeLoad::fem = fem;
		RD = FeModel::RD;

		spLoad = std::vector<double>(fem->nEq);
		dpLoad = std::vector<double>(fem->nEq);
		dhLoad = std::vector<double>(fem->nEq);
		sDispl = std::vector<double>(fem->nEq);
		dDispl = std::vector<double>(fem->nEq);
		RHS = std::vector<double>(fem->nEq);

		if (fem->thermalLoading)
		{
			dtemp = std::vector<double>(fem->nNod);
		}
	}

	bool FeLoad::readData()
	{

		return readDataFile(RD, true);
	}

	bool FeLoad::readDataFile(FeScanner *es, bool newLoad)
	{
		if (newLoad)
		{
			scaleLoad = 0;
			nodForces = std::list();
			itnf = nodForces.listIterator(0);
			surForces = std::list();
			itsf = surForces.listIterator(0);
			if (fem->thermalLoading)
			{
				for (int i = 0; i < dtemp.size(); i++)
				{
					dtemp[i] = 0.0;
				}
			}
			for (int i = 0; i < dDispl.size(); i++)
			{
				dDispl[i] = 0;
			}
		}

		if (!es->hasNext())
		{
			return false; // No load data
		}

		vars name = nullptr;
		std::wstring s;

		while (es->hasNext())
		{
			std::wstring varName = es->next();
			std::wstring varNameLower = StringHelper::toLower(varName);
			if (varName == L"#")
			{
				es->nextLine();
				continue;
			}
			try
			{
				name = varsHelper::enumFromString(varNameLower);
			}
			catch (const std::runtime_error &e)
			{
				UTIL::errorMsg(L"Variable name is not found: " + varName);
			}

			switch (name)
			{

			case model::FeLoadData::vars::loadstep:
				loadStepName = es->next();
				break;

			case model::FeLoadData::vars::scaleload:
				scaleLoad = es->readDouble();
				break;

			case model::FeLoadData::vars::residtolerance:
				residTolerance = es->readDouble();
				break;

			case model::FeLoadData::vars::maxiternumber:
				maxIterNumber = es->readInt();
				break;

			case model::FeLoadData::vars::nodforce:
				readNodalForces(es);
				break;

			case model::FeLoadData::vars::surforce:
				readSurForces(es);
				break;

			case model::FeLoadData::vars::boxsurforce:
				createBoxSurForces(es);
				break;

			case model::FeLoadData::vars::nodtemp:
				dtemp = std::vector<double>(fem->nNod);
				for (int i = 0; i < fem->nNod; i++)
				{
					dtemp[i] = es->readDouble();
				}
				break;

			case model::FeLoadData::vars::includefile:
			{
				s = StringHelper::toLower(es->next());
				FeScanner *R = new FeScanner(s);
				readDataFile(R, false);

//JAVA TO C++ CONVERTER TODO TASK: A 'delete R' statement was not added since R was passed to a method or constructor. Handle memory management manually.
				break;

//JAVA TO C++ CONVERTER TODO TASK: A 'delete R' statement was not added since R was passed to a method or constructor. Handle memory management manually.
			}
			case model::FeLoadData::vars::end:
				return true;
			}
			es++;
		}
		return true;
	}

	void FeLoad::readNodalForces(FeScanner *es)
	{

		std::wstring s = StringHelper::toLower(es->next());
		int idf = UTIL::direction(s);
		if (idf == -1)
		{
			UTIL::errorMsg(L"nodForce" + L" direction should be x/y/z. Specified:" + s);
		}

		if (!es->hasNextDouble())
		{
			UTIL::errorMsg(L"nodForce value is not a double: " + es->next());
		}
		double vd = es->nextDouble();

		itnf = es->readNumberList(itnf, idf, fem->nDim, vd);
	}

	void FeLoad::readSurForces(FeScanner *es)
	{

		std::wstring s = StringHelper::toLower(es->next());
		int dir = UTIL::direction(s);
		if (dir == -1)
		{
			UTIL::errorMsg(L"surForce" + L" direction should be x/y/z/n. Specified:" + s);
		}
		int iel = es->readInt();
		int nFaceNodes = es->readInt();
		for (int i = 0; i < nFaceNodes; i++)
		{
			iw[i] = es->readInt();
		}
		for (int i = 0; i < nFaceNodes; i++)
		{
			dw[i] = es->readDouble();
		}
		ElemFaceLoad tempVar(iel - 1,nFaceNodes,dir,iw,dw);
		itsf->add(&tempVar);
	}

	void FeLoad::createBoxSurForces(FeScanner *es)
	{
		std::vector<std::vector<int>> faces;
		std::wstring s = StringHelper::toLower(es->next());
		int dir = UTIL::direction(s);
		if (dir == -1)
		{
			UTIL::errorMsg(L"boxSurForce" + L" direction should be x/y/z/n. Specified:" + s);
		}

		if (!es->hasNextDouble())
		{
			UTIL::errorMsg(L"boxSurForce value is not a double: " + es->next());
		}
		double force = es->nextDouble();

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < fem->nDim; j++)
			{
				box[i][j] = es->readDouble();
			}
		}

		for (int iel = 0; iel < fem->nEl; iel++)
		{
			Element *el = fem->elems[iel];
			faces = el->getElemFaces();
			for (auto face : faces)
			{
				int nNodes = face.size();
				for (int inod = 0; inod < nNodes; inod++)
				{
					iw[inod] = 0;
				}
				for (int inod = 0; inod < nNodes; inod++)
				{
					int iGl = el->ind[face[inod]];
					if (iGl > 0)
					{
						for (int j = 0; j < fem->nDim; j++)
						{
							double x = fem->getNodeCoord(iGl - 1,j);
							if (x < box[0][j] || x > box[1][j])
							{
								goto FACEContinue;
							}
						}
						iw[inod] = iGl;
					}
				}
				ElemFaceLoad tempVar(iel,nNodes,dir,iw,force);
				itsf->add(&tempVar);
				FACEContinue:;
			}
			FACEBreak:;
		}
	}

	void FeLoad::assembleRHS()
	{

		if (scaleLoad != 0.0)
		{
			for (int i = 0; i < fem->nEq; i++)
			{
				dpLoad[i] *= scaleLoad;
				dhLoad[i] *= scaleLoad;
				RHS [i] = dpLoad[i] + dhLoad[i];
			}
			return;
		}
		for (int i = 0; i < fem->nEq; i++)
		{
			dpLoad[i] = 0.0;
			dhLoad[i] = 0.0;
		}

		// Nodal forces specified directly
		itnf = nodForces.listIterator(0);
		Dof *d;
		while (itnf->hasNext())
		{
			d = static_cast<Dof*>(itnf->next());
			dpLoad[d->dofNum - 1] = d->value;
			itnf++;
		}

		// Surface load at element faces
		itsf = surForces.listIterator(0);
		ElemFaceLoad *efl;
		Element *elm;
		while (itsf->hasNext())
		{
			efl = static_cast<ElemFaceLoad*>(itsf->next());
			elm = fem->elems[efl->iel];
			elm->setElemXy();
			if (elm->equivFaceLoad(efl) == -1)
			{
				UTIL::errorMsg(L"surForce" + L" does not match any face of element: " + std::to_wstring(efl->iel));
			}
			elm->assembleElemVector(Element::evec,dpLoad);
			itsf++;
		}

		// Temperature field
		if (fem->thermalLoading)
		{
			for (int iel = 0; iel < fem->nEl; iel++)
			{
				elm = fem->elems[iel];
				elm->setElemXyT();
				elm->thermalVector();
				elm->assembleElemVector(Element::evec,dhLoad);
			}
		}

		// Right-hand side = actual load + fictitious load
		for (int i = 0; i < fem->nEq; i++)
		{
			RHS[i] = dpLoad[i] + dhLoad[i];
		}

		// Displacement boundary conditions for right-hand side
		ListIterator itdbc = fem->defDs.listIterator(0);
		while (itnf->hasNext())
		{
			d = static_cast<Dof*>(itdbc->next());
			RHS[d->dofNum - 1] = FE::bigValue * d->value;
			itnf++;
		}
	}
}
