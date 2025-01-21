#include "FeModelData.h"
//#include "../util/FeScanner.h"
#include "../elem/Element.h"

namespace model
{
	using namespace elem;
	using namespace util;
    FeScanner *FeModelData::RD;
    //java::io::PrintWriter *FeModelData::PR;
    std::wstring FeModelData::varName;
    FeModelData::StrStates FeModelData::stressState = StrStates::threed;

	void FeModelData::newCoordArray()
	{
		xyz = std::vector<double>(nNod * nDim);
	}

	void FeModelData::setNodeCoords(int node, std::vector<double> &xyzn)
	{
		for (int i = 0; i < nDim; i++)
		{
			xyz[node * nDim + i] = xyzn[i];
		}
	}

	void FeModelData::setNodeCoord(int node, int i, double v)
	{
		xyz[node * nDim + i] = v;
	}

	std::vector<double> FeModelData::getNodeCoords(int node)
	{
		std::vector<double> nodeCoord(nDim);
		for (int i = 0; i < nDim; i++)
		{
			nodeCoord[i] = xyz[node * nDim + i];
		}
		return nodeCoord;
	}

	double FeModelData::getNodeCoord(int node, int i)
	{
		return xyz[node * nDim + i];
	}
}
