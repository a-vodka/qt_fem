#include "Material.h"
#include "ElasticMaterial.h"
#include "ElasticPlasticMaterial.h"
#include "../elem/Element.h"

namespace material
{
	using Element = elem::Element;

    Material *Material::newMaterial(const wchar_t* matPhysLaw, const wchar_t* stressState)
    {
        if (wcscmp(matPhysLaw, L"elastic") == 0)
		{
			  return new ElasticMaterial(stressState);
		}
		else
        {
              //return new ElasticPlasticMaterial(std::wstring(stressState));
              return nullptr;
		}
	}

	void Material::strainToStress(Element *elm, int ip)
	{
	}

	void Material::setElasticProp(double e, double nu, double alpha)
	{
		this->e = e;
		this->nu = nu;
		this->alpha = alpha;
	}

	void Material::setPlasticProp(double sY, double km, double mm)
	{
		this->sY = sY;
		this->km = km;
		this->mm = mm;
	}

	double Material::getLambda()
	{
		return (stressState == L"plstress") ? e * nu / ((1 + nu) * (1 - nu)) : e * nu / ((1 + nu) * (1 - 2 * nu));
	}

	double Material::getMu()
	{
		return 0.5 * e / (1 + nu);
	}

	double Material::getNu()
	{
		return nu;
	}

	double Material::getAlpha()
	{
		return alpha;
	}

	void Material::elasticityMatrix(std::vector<std::vector<double>> &emat)
	{
	}
}
