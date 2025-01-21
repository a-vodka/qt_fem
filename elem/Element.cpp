#include "Element.h"
//#include "../model/FeModel.h"
#include "../model/FeLoad.h"
#include "../material/Material.h"
#include "StressContainer.h"
#include "../util/UTIL.h"
#include "../fea/FE.h"
#include "../model/ElemFaceLoad.h"
#include "elem/ElementQuad3D.h"
#include "model/FeModel.h"

namespace elem
{
	using namespace model;
	using namespace material;
	using FE = fea::FE;
	using UTIL = util::UTIL;
    FeModel *Element::fem;
    FeLoadData *Element::load;
    Material *Element::mat;
    std::vector<std::vector<double>> Element::kmat = RectangularVectors::RectangularDoubleVector(60, 60);
    std::vector<double> Element::evec(60);
    std::vector<std::vector<double>> Element::xy = RectangularVectors::RectangularDoubleVector(20, 3);
    std::vector<double> Element::dtn(20);
    std::vector<double> Element::dstrain(6);

	Element *Element::newElement(const std::wstring &name)
	{

/*        elements* el = nullptr;
		try
        {
            el = elementsHelper::enumFromString(name);
		}
		catch (const std::runtime_error &e)
		{
			UTIL::errorMsg(L"Incorrect element type: " + name);
        }*/

        return new ElementQuad3D();
	}

	Element::Element(const std::wstring &name, int nind, int nstress)
	{
		this->name = name;
		ind = std::vector<int>(nind);

        str = std::vector<StressContainer*>(nstress);
        for (int ip = 0; ip < nstress; ip++)
        {
            str[ip] = new StressContainer(fem->nDim);
        }

	}

	void Element::stiffnessMatrix()
	{
	}

	void Element::thermalVector()
	{
	}

	int Element::equivFaceLoad(ElemFaceLoad *surLd)
	{
		return -1;
	}

	void Element::equivStressVector()
	{
	}

	std::vector<std::vector<int>> Element::getElemFaces()
	{
		return std::vector<std::vector<int>>
		{
			std::vector<int> {0},
			std::vector<int> {0}
		};
	}

	std::vector<double> Element::getStrainsAtIntPoint(int intPoint)
	{
		return std::vector<double> {0, 0};
	}

	double Element::getTemperatureAtIntPoint(int intPoint)
	{
		return 0.0;
	}

	void Element::extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn)
	{
	}

	void Element::setElemConnectivities(std::vector<int> &indel, int nind)
	{
		std::copy_n(indel.begin(), nind, ind.begin());
	}

	void Element::setElemConnectivities(std::vector<int> &indel)
	{
		std::copy_n(indel.begin(), indel.size(), ind.begin());
	}

	void Element::setElemMaterial(const std::wstring &mat)
	{
		matName = mat;
	}

	void Element::setElemXy()
	{
        for (size_t i = 0; i < ind.size(); i++)
		{
			int indw = ind[i] - 1;
			if (indw >= 0)
			{
				xy[i] = fem->getNodeCoords(indw);
			}
		}
	}

	void Element::setElemXyT()
	{
        for (size_t i = 0; i < ind.size(); i++)
		{
			int indw = ind[i] - 1;
			if (indw >= 0)
			{
				if (fem->thermalLoading)
				{
					dtn[i] = FeLoad::dtemp[indw];
				}
				xy[i] = fem->getNodeCoords(indw);
			}
		}
	}

	void Element::assembleElemVector(std::vector<double> &elVector, std::vector<double> &glVector)
	{
        for (size_t i = 0; i < ind.size(); i++)
		{
			int indw = ind[i] - 1;
			if (indw >= 0)
			{
				int adr = indw * fem->nDim;
				for (int j = 0; j < fem->nDim; j++)
				{
					glVector[adr + j] += elVector[i * fem->nDim + j];
				}
			}
		}
	}

	void Element::disAssembleElemVector(std::vector<double> &glVector)
	{
        for (size_t i = 0; i < ind.size(); i++)
		{
			int indw = ind[i] - 1;
			if (indw >= 0)
			{
				int adr = indw * fem->nDim;
				for (int j = 0; j < fem->nDim; j++)
				{
					evec[i * fem->nDim + j] = glVector[adr + j];
				}
			}
		}
	}

	std::vector<int> Element::getElemConnectivities()
	{
		std::vector<int> indE(ind.size());
		std::copy_n(ind.begin(), ind.size(), indE.begin());
		return indE;
	}

	void Element::accumulateStress()
	{
        for (size_t ip = 0; ip < str.size(); ip++)
		{
			for (int i = 0; i < 2 * fem->nDim; i++)
			{
				str[ip]->sStress[i] += str[ip]->dStress[i];
			}
		}
            for (size_t ip = 0; ip < str.size(); ip++)
			{
				str[ip]->sEpi += str[ip]->dEpi;
			}
	}
}
