#include "FeModel.h"
//#include "../util/FeScanner.h"
#include "../material/Material.h"
#include "../util/UTIL.h"
#include "../solver/Solver.h"
#include "../elem/Element.h"
#include "../model/Dof.h"

namespace model
{
	using namespace elem;
	using namespace material;
	using namespace util;
	using namespace solver;
    std::vector<int> FeModel::elCon(20);
    std::vector<std::vector<double>> FeModel::box = RectangularVectors::RectangularDoubleVector(2, 3);
    /*
	FeModel::FeModel(FeScanner *RD, PrintWriter *PR)
	{
		FeModelData::RD = RD;
		FeModelData::PR = PR;
	}
    */
    /*
	void FeModel::readData()
	{
		readDataFile(RD);
	}

	void FeModel::readDataFile(FeScanner *es)
	{
/*
		vars name = nullptr;
		std::wstring s;
		Material *mat;
		it = defDs.listIterator(0);

		while (es->hasNext())
		{
			varName = es->next();
			std::wstring varname = StringHelper::toLower(varName);
			if (varName == L"#")
			{
				es->nextLine();
				continue;
			}
			try
			{
				name = varsHelper::enumFromString(varname);
			}
			catch (const std::runtime_error &e)
			{
				UTIL::errorMsg(L"Variable name is not found: " + varName);
			}

			switch (name)
			{

			case model::FeModelData::vars::nel:
				nEl = es->readInt();
				break;

			case model::FeModelData::vars::nnod:
				nNod = es->readInt();
				break;

			case model::FeModelData::vars::ndim:
				nDim = es->readInt();
				nDf = nDim;
				break;

			case model::FeModelData::vars::stressstate:
				s = StringHelper::toLower(es->next());
				try
				{
					stressState = StrStatesHelper::enumFromString(s);
				}
				catch (const std::runtime_error &e)
				{
					UTIL::errorMsg(L"stressState has forbidden value: " + s);
				}
				if (stressState != StrStates::threed)
				{
					  nDim = nDf = 2;
				}
				else
				{
					nDim = nDf = 3;
				}
				break;

			case model::FeModelData::vars::physlaw:
				s = StringHelper::toLower(es->next());
				try
				{
					physLaw = PhysLawsHelper::enumFromString(s);
				}
				catch (const std::runtime_error &e)
				{
					UTIL::errorMsg(L"physLaw has forbidden value: " + s);
				}
				break;

			case model::FeModelData::vars::solver:
				s = StringHelper::toLower(es->next());
				try
				{
					Solver::solver = Solver::SolversHelper::enumFromString(s);
				}
				catch (const std::runtime_error &e)
				{
					UTIL::errorMsg(L"solver has forbidden value: " + s);
				}
				break;

			case model::FeModelData::vars::elcon:
				readElemData(es);
				break;

			case model::FeModelData::vars::nodcoord:
				if (nNod == 0 || nDim == 0)
				{
					UTIL::errorMsg(L"nNod and nDim should be" + L" specified before nodCoord");
				}
				nEq = nNod * nDim;
				// Nodal coordinates
				newCoordArray();
				for (int i = 0; i < nNod; i++)
				{
					for (int j = 0; j < nDim; j++)
					{
						setNodeCoord(i, j, es->readDouble());
					}
				}
				break;

			case model::FeModelData::vars::material:
			{
				std::wstring matname = es->next();
				mat = Material::newMaterial(model::FeModelData::PhysLawsHelper::enumName(physLaw), model::FeModelData::StrStatesHelper::enumName(stressState));
				double e = es->readDouble();
				double nu = es->readDouble();
				double alpha = es->readDouble();
				mat->setElasticProp(e, nu, alpha);
				if (physLaw == PhysLaws::elplastic)
				{
					double sY = es->readDouble();
					double km = es->readDouble();
					double mm = es->readDouble();
					mat->setPlasticProp(sY, km, mm);
				}
				materials.emplace(matname, mat);
				break;

			}
			case model::FeModelData::vars::constrdispl:
				readConstrDisplacements(es);
				break;

			case model::FeModelData::vars::boxconstrdispl:
				createBoxConstrDisplacements(es);
				break;

			case model::FeModelData::vars::thermalloading:
				s = es->next();
				if (StringHelper::toLower(s).equals(L"y"))
				{
					thermalLoading = true;
				}
				else if (StringHelper::toLower(s).equals(L"n"))
				{
					thermalLoading = false;
				}
				else
				{
					UTIL::errorMsg(L"thermalLoading should be" + L" y/n. Specified: " + s);
				}
				break;

			case model::FeModelData::vars::includefile:
			{
				s = StringHelper::toLower(es->next());
				FeScanner *R = new FeScanner(s);
				readDataFile(R);

//JAVA TO C++ CONVERTER TODO TASK: A 'delete R' statement was not added since R was passed to a method or constructor. Handle memory management manually.
				break;

//JAVA TO C++ CONVERTER TODO TASK: A 'delete R' statement was not added since R was passed to a method or constructor. Handle memory management manually.
			}
			case model::FeModelData::vars::end:
				return;
			}
			es++;
		}
*/
	}
/*
	void FeModel::readElemData(FeScanner *es)
	{

		if (nEl == 0)
		{
			UTIL::errorMsg(L"nEl should be defined before elCon");
		}
		elems = std::vector<Element*>(nEl);

		for (int iel = 0; iel < nEl; iel++)
		{
			// Element type
			std::wstring s = StringHelper::toLower(es->next());
			elems[iel] = Element::newElement(s);
			// Element material
			std::wstring elMat = es->next();
			elems[iel]->setElemMaterial(elMat);
			// Element connectivities
			int nind = elems[iel]->ind->size();
			for (int i = 0; i < nind; i++)
			{
				elCon[i] = es->readInt();
			}
			elems[iel]->setElemConnectivities(elCon,nind);
		}
	}

	void FeModel::readConstrDisplacements(FeScanner *es)
	{
		std::wstring s = StringHelper::toLower(es->next());
		int idf = UTIL::direction(s);
		if (idf == -1)
		{
			UTIL::errorMsg(L"constrDispl direction" + L" should be x/y/z. Specified:" + s);
		}
		if (!es->hasNextDouble())
		{
			UTIL::errorMsg(L"constrDispl value is not a double: " + es->next());
		}
		double vd = es->nextDouble();
		it = es->readNumberList(it, idf, nDim, vd);
	}

	void FeModel::createBoxConstrDisplacements(FeScanner *es)
	{
		std::wstring s = StringHelper::toLower(es->next());
		int idf = UTIL::direction(s);
		if (idf == -1)
		{
			UTIL::errorMsg(L"boxConstrDispl direction should be" + L" x/y/z. Specified:" + s);
		}
		if (!es->hasNextDouble())
		{
			UTIL::errorMsg(L"boxConstrDispl value is not" + L" a double: " + es->next());
		}
		double vd = es->nextDouble();
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < nDim; j++)
			{
				box[i][j] = es->readDouble();
			}
		}
		for (int i = 0; i < nNod; i++)
		{
			for (int j = 0; j < nDim; j++)
			{
				double x = getNodeCoord(i,j);
				if (x < box[0][j] || x> box[1][j])
				{
					goto nodeContinue;
				}
			}
			Dof tempVar(nDim * i + idf, vd);
			it->add(&tempVar);
			nodeContinue:;
		}
		nodeBreak:;

	}

}
*/
