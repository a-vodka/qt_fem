#include "ElasticPlasticMaterial.h"
#include "../util/UTIL.h"
#include "../fea/FE.h"
#include "../elem/Element.h"
#include "elem/StressContainer.h"

namespace material
{
	using Element = elem::Element;
	using FE = fea::FE;
	using UTIL = util::UTIL;
    double sig0[6];
    double sig[6];
    double a[6];
    double Emat[6][6];
    double Ea[6];
    double depsp[6];
    double ElasticPlasticMaterial::epi = 0;
    double ElasticPlasticMaterial::G = 0;

	ElasticPlasticMaterial::ElasticPlasticMaterial(const std::wstring &stressState) : ElasticMaterial(stressState)
	{

		if (stressState == L"plstress")
		{
            UTIL::errorMsg(L"Elastic-plastic material is not implemented for plane stress");
		}
		if (!FE::epIntegrationTANGENT)
		{
			midpoint = new Midpoint(this);
		}
	}

	void ElasticPlasticMaterial::strainToStress(Element *elm, int ip)
	{

		elasticityMatrix(Emat);
		G = getMu();

		// Elastic stress increment dsig due to deps
        this->strainToStress(elm, ip);

		for (int i = 0; i < lv; i++)
		{
			sig0[i] = elm->str[ip]->sStress[i];
			sig[i] = sig0[i] + dsig[i];
		}
		// Equivalent plastic strain
		epi = elm->str[ip]->sEpi;
		double epi0 = epi;
		double f1 = yieldFunction(sig, epi);

		// Elastic point
		if (f1 < 0.0)
		{
			return;
		}

		// Elastic-plastic point
		double f0 = yieldFunction(sig0, epi);
		double r;
		if (f0 < 0.0)
		{
			r = -f0 / (f1 - f0);
			for (int i = 0; i < lv; i++)
			{
				sig[i] = sig0[i] + dsig[i] * r;
			}
			double f = yieldFunction(sig, epi);
			derivYieldFunc(sig, a);
			double c1 = 0.0;
			for (int i = 0; i < lv; i++)
			{
				c1 += a[i] * dsig[i];
			}
			r = r - f / c1;
		}
		else
		{
			r = 0.0;
		}

		// Number of subincrements ( = 1 for midpoint method)
		int nsub = (FE::epIntegrationTANGENT) ? static_cast<int>(f1 / (beta * sY)) + 1 : 1;

		for (int i = 0; i < lv; i++)
		{
			sig[i] = sig0[i] + dsig[i] * r;
			dsig[i] = (1.0 - r) * dsig[i] / nsub;
			deps[i] = (1.0 - r) * deps[i] / nsub;
		}
		// Subincrement loop: tangent or midpoint method
		for (int isub = 0; isub < nsub; isub++)
		{
			if (FE::epIntegrationTANGENT)
			{
				 tangentStressIncrement();
			}
			else
			{
				midpoint->stressIncrement();
			}
		}
		for (int i = 0; i < lv; i++)
		{
			elm->str[ip]->dStress[i] = sig[i] - sig0[i];
		}
		elm->str[ip]->dEpi = epi - epi0;
	}

	void ElasticPlasticMaterial::tangentStressIncrement()
	{

		double H = slopeH(epi);
		derivYieldFunc(sig, a);
		double dlambda = 0.0;
		for (int i = 0; i < lv; i++)
		{
			double s = 0.0;
			for (int j = 0; j < lv; j++)
			{
				s += Emat[i][j] * a[j];
			}
			Ea[i] = s;
			dlambda += a[i] * dsig[i];
		}
		dlambda /= (H + 3 * G);
		if (dlambda < 0.0)
		{
			dlambda = 0.0;
		}
		for (int i = 0; i < lv; i++)
		{
			sig[i] += dsig[i] - dlambda * Ea[i];
			depsp[i] = dlambda * a[i];
		}
		epi += eqPlastStrain(depsp);

		// Stress correction
		double f = yieldFunction(sig, epi);
		double c1 = 0.0;
		for (int i = 0; i < lv; i++)
		{
			c1 += a[i] * a[i];
		}
		for (int i = 0; i < lv; i++)
		{
			sig[i] -= a[i] * f / c1;
		}
	}

	double ElasticPlasticMaterial::yieldFunction(std::vector<double> &s, double epi)
	{
		double sm, seq;
		if (stressState == L"threed")
		{
			sm = (s[0] + s[1] + s[2]) / 3;
			seq = std::sqrt(3.0 * (0.5 * ((s[0] - sm) * (s[0] - sm) + (s[1] - sm) * (s[1] - sm) + (s[2] - sm) * (s[2] - sm)) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
		}
		else
		{
			sm = (s[0] + s[1] + s[3]) / 3;
			seq = std::sqrt(3.0 * (0.5 * ((s[0] - sm) * (s[0] - sm) + (s[1] - sm) * (s[1] - sm) + (s[3] - sm) * (s[3] - sm)) + s[2] * s[2]));
		}
		return seq - yieldRadius(epi);
	}

	double ElasticPlasticMaterial::yieldRadius(double ep)
	{
		if (ep <= 0.0)
		{
			return sY;
		}
		else
		{
			return (km * std::pow(ep, mm) + sY);
		}
	}

	void ElasticPlasticMaterial::derivYieldFunc(std::vector<double> &s, std::vector<double> &a)
	{
		double sm, seq;

		if (stressState == L"threed")
		{
			sm = (s[0] + s[1] + s[2]) / 3;
			a[0] = s[0] - sm;
			a[1] = s[1] - sm;
			a[2] = s[2] - sm;
			a[3] = 2 * s[3];
			a[4] = 2 * s[4];
			a[5] = 2 * s[5];
			seq = std::sqrt(0.5 * (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) + s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
		}
		else
		{
			sm = (s[0] + s[1] + s[3]) / 3;
			a[0] = s[0] - sm;
			a[1] = s[1] - sm;
			a[2] = 2 * s[2];
			a[3] = s[3] - sm;
			seq = std::sqrt(0.5 * (a[0] * a[0] + a[1] * a[1] + a[3] * a[3]) + s[2] * s[2]);
		}

		for (int i = 0; i < lv; i++)
		{
			a[i] = 0.5 * std::sqrt(3.0) / seq * a[i];
		}
	}

	double ElasticPlasticMaterial::slopeH(double epi)
	{
		if (km == 0.0)
		{
			return 0.0;
		}
		else if (mm == 1.0)
		{
			return km;
		}
		else
		{
			return (epi == 0.0) ? 0.0 : km * mm * std::pow(epi, mm - 1.0);
		}
	}

	double ElasticPlasticMaterial::eqPlastStrain(std::vector<double> &dp)
	{
		if (stressState == L"threed")
		{
			return std::sqrt((2 * (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]) + dp[3] * dp[3] + dp[4] * dp[4] + dp[5] * dp[5]) / 3.0);
		}
		else
		{
			if (stressState == L"plstress")
			{
				dp[3] = -(dp[0] + dp[1]);
			}
			return std::sqrt((2 * (dp[0] * dp[0] + dp[1] * dp[1] + dp[3] * dp[3]) + dp[2] * dp[2]) / 3.0);
		}
	}

	ElasticPlasticMaterial::Midpoint::Midpoint(ElasticPlasticMaterial *outerInstance) : outerInstance(outerInstance)
	{
	}

	void ElasticPlasticMaterial::Midpoint::stressIncrement()
	{
		double SQ32 = std::sqrt(1.5);
		double SQ23 = std::sqrt(2.0 / 3.0);
		double tolerance = 1.e-5;
		// Transform strains to tensor components
		if (lv == 4)
		{
			deps[2] *= 0.5;
		}
		else
		{
			for (int i = 3; i < 6; i++)
			{
				deps[i] *= 0.5;
			}
		}
		double depsm = deviator(deps, ed);
		double sigm = deviator(sig, sd);
		double sigeq0 = SQ32 * norm(sd);
		for (int i = 0; i < lv; i++)
		{
			sdtr[i] = sd[i] + 2 * G * ed[i];
			a0[i] = 1.5 * sd[i] / sigeq0;
		}
		double lambda = SQ23 * norm(ed);
		// Find lambda by Newton-Raphson iteration
		double epi1, sigeq;
		for (; ;)
		{
			for (int i = 0; i < lv; i++)
			{
				sal[i] = sdtr[i] - 2 * G * lambda * (1 - alpha) * a0[i];
			}
			double salmod = norm(sal);
			for (int i = 0; i < lv; i++)
			{
				salbar[i] = sal[i] / salmod;
				b[i] = (1 - alpha) * a0[i] + alpha * SQ32 * salbar[i];
			}
			double bmod = norm(b);
			epi1 = epi + lambda * SQ23 * bmod;
			sigeq = outerInstance->yieldRadius(epi1);
			double phi = SQ32 * salmod - 3 * G * alpha * lambda - sigeq;
			if (std::abs(phi) < tolerance * outerInstance->sY)
			{
				break;
			}
			double phiPrime = -2.0 * SQ32 * G * (1 - alpha) * dyadicProduct(a0, salbar) - 3.0 * G * alpha - outerInstance->slopeH(epi1) * SQ23 * bmod;
			double lambda1 = lambda - phi / phiPrime;
			lambda = (lambda1 <= 0.0) ? 0.5 * lambda:lambda1;
		}
		epi = epi1;
		for (int i = 0; i < lv; i++)
		{
			sd[i] = sal[i] / (1 + 3 * G * alpha * lambda / sigeq);
		}
		double dsigm = depsm * outerInstance->e / (1.0 - 2.0 * outerInstance->nu);
		stressFromDeviator(sd, sigm + dsigm, sig);
	}

    double ElasticPlasticMaterial::Midpoint::deviator(double s[], double d[])
	{
		double sm;
		if (lv == 4)
		{
			sm = (s[0] + s[1] + s[3]) / 3;
			d[0] = s[0] - sm;
			d[1] = s[1] - sm;
			d[3] = s[3] - sm;
			d[2] = s[2];
		}
		else
		{
			sm = (s[0] + s[1] + s[2]) / 3;
			d[0] = s[0] - sm;
			d[1] = s[1] - sm;
			d[2] = s[2] - sm;
			d[3] = s[3];
			d[4] = s[4];
			d[5] = s[5];
		}
		return (sm);
	}
    double ElasticPlasticMaterial::Midpoint::deviator(std::vector<double> &s, double d[])
    {
        double sm;
        if (lv == 4)
        {
            sm = (s[0] + s[1] + s[3]) / 3;
            d[0] = s[0] - sm;
            d[1] = s[1] - sm;
            d[3] = s[3] - sm;
            d[2] = s[2];
        }
        else
        {
            sm = (s[0] + s[1] + s[2]) / 3;
            d[0] = s[0] - sm;
            d[1] = s[1] - sm;
            d[2] = s[2] - sm;
            d[3] = s[3];
            d[4] = s[4];
            d[5] = s[5];
        }
        return (sm);
    }

	void ElasticPlasticMaterial::Midpoint::stressFromDeviator(std::vector<double> &d, double sm, std::vector<double> &s)
	{
		if (lv == 4)
		{
			s[0] = d[0] + sm;
			s[1] = d[1] + sm;
			s[3] = d[3] + sm;
			s[2] = d[2];
		}
		else
		{
			s[0] = d[0] + sm;
			s[1] = d[1] + sm;
			s[2] = d[2] + sm;
			s[3] = d[3];
			s[4] = d[4];
			s[5] = d[5];
		}
	}

    double ElasticPlasticMaterial::Midpoint::norm(double s[])
	{
		if (lv == 4)
		{
			return (std::sqrt(s[0] * s[0] + s[1] * s[1] + s[3] * s[3] + 2 * s[2] * s[2]));
		}
		else
		{
			return (std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])));
		}
	}

	double ElasticPlasticMaterial::Midpoint::dyadicProduct(std::vector<double> &a, std::vector<double> &b)
	{
		if (lv == 4)
		{
			return (a[0] * b[0] + a[1] * b[1] + a[3] * b[3] + 2 * a[2] * b[2]);
		}
		else
		{
			return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]));
		}
	}
}
