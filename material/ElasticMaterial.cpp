#include "ElasticMaterial.h"
#include "../elem/Element.h"
#include "elem/StressContainer.h"

namespace material
{
    using Element = elem::Element;
    std::vector<double> ElasticMaterial::deps;
    std::vector<double> ElasticMaterial::dsig;
    //double deps[6];
    //double dsig[6];
    int ElasticMaterial::lv = 0;

	ElasticMaterial::ElasticMaterial(const std::wstring &stressState)
	{

		this->stressState = stressState;
		lv = (stressState == L"threed")? 6:4;
        deps.resize(lv);
        dsig.resize(lv);
	}

	void ElasticMaterial::strainToStress(Element *elm, int ip)
	{

		deps = elm->getStrainsAtIntPoint(ip);
		double temp = elm->getTemperatureAtIntPoint(ip);

		double mu = 0.5 * e / (1 + nu);
		double lambda = (stressState == L"plstress") ? e * nu / (1 - nu * nu) : e * nu / ((1 + nu) * (1 - 2 * nu));
		double beta = lambda + 2.0 * mu;
		double at = alpha * temp;

		if (stressState == L"threed")
		{
			deps[0] -= at;
			deps[1] -= at;
			deps[2] -= at;
			dsig[0] = beta * deps[0] + lambda * (deps[1] + deps[2]);
			dsig[1] = beta * deps[1] + lambda * (deps[0] + deps[2]);
			dsig[2] = beta * deps[2] + lambda * (deps[0] + deps[1]);
			dsig[3] = mu * deps[3];
			dsig[4] = mu * deps[4];
			dsig[5] = mu * deps[5];
		}
		else
		{
			deps[0] -= at;
			deps[1] -= at;
			if (stressState != L"plstress")
			{
				deps[3] -= at;
			}
			dsig[0] = beta * deps[0] + lambda * (deps[1] + deps[3]);
			dsig[1] = beta * deps[1] + lambda * (deps[0] + deps[3]);
			dsig[2] = mu * deps[2];
			dsig[3] = 0.0;
			if (stressState == L"plstrain")
			{
			  dsig[3] = nu * (dsig[0] + dsig[1]) - e * at;
			}
			if (stressState == L"axisym")
			{
			  dsig[3] = beta * deps[3] + lambda * (deps[0] + deps[1]);
			}
		}
		for (int i = 0; i < lv; i++)
		{
            elm->str[ip]->dStress[i] = dsig[i];
		}
	}

	void ElasticMaterial::elasticityMatrix(std::vector<std::vector<double>> &emat)
	{
		if (stressState == L"threed")
		{
			elasticityMatrix3D(emat);
		}
		else
		{
			elasticityMatrix2D(emat);
		}
	}

	void ElasticMaterial::elasticityMatrix3D(std::vector<std::vector<double>> &emat)
	{

		double mu = getMu();
		double lambda = getLambda();
		double beta = lambda + 2 * mu;

		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				emat[i][j] = 0;
			}
		}

		emat[0][0] = emat[1][1] = emat[2][2] = beta;
		emat[0][1] = emat[1][0] = emat[0][2] = emat[2][0] = emat[1][2] = emat[2][1] = lambda;
		emat[3][3] = emat[4][4] = emat[5][5] = mu;
	}

	void ElasticMaterial::elasticityMatrix2D(std::vector<std::vector<double>> &emat)
	{

		double mu = getMu();
		double lambda = getLambda();
		double beta = lambda + 2 * mu;

		emat[0][0] = emat[3][3] = emat[1][1] = beta;
		emat[0][1] = emat[1][0] = emat[0][3] = emat[3][0] = emat[1][3] = emat[3][1] = lambda;
		emat[2][2] = mu;
		emat[0][2] = emat[2][0] = emat[1][2] = emat[2][1] = emat[2][3] = emat[3][2] = 0.0;
	}
}
