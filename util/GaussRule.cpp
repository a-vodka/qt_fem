#include "GaussRule.h"
#include "UTIL.h"
#include <cmath>
namespace util
{

const std::vector<std::vector<double>> GaussRule::X =
{
	{0.0},
	{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)},
	{-std::sqrt(0.6), 0.0, std::sqrt(0.6)}
};
const std::vector<std::vector<double>> GaussRule::W =
{
	{2.0},
	{1.0, 1.0},
	{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}
};
const std::vector<double> GaussRule::X14 = {-a, a, -a, -a, a, -a, a, a, -b, b, 0, 0, 0, 0};
const std::vector<double> GaussRule::Y14 = {-a, -a, a, -a, a, a, -a, a, 0, 0, -b, b, 0, 0};
const std::vector<double> GaussRule::Z14 = {-a, -a, -a, a, -a, a, a, a, 0, 0, 0, 0, -b, b};

	GaussRule::GaussRule(int nGauss, int nDim)
	{

		if (!((nGauss >= 1 && nGauss <= 3) || nGauss == 14))
		{
			UTIL::errorMsg(L"nGauss has forbidden value: " + std::to_wstring(nGauss));
		}
		if (!(nDim >= 1 && nDim <= 3))
		{
			UTIL::errorMsg(L"GaussRule: nDim has forbidden value: " + std::to_wstring(nDim));
		}

		if (nGauss == 14)
		{
			nIntPoints = 14;
		}
		else
		{
			nIntPoints = 1;
			for (int i = 0; i < nDim; i++)
			{
				nIntPoints *= nGauss;
			}
		}

		xii = std::vector<double>(nIntPoints);
		wi = std::vector<double>(nIntPoints);
		if (nDim > 1)
		{
			eti = std::vector<double>(nIntPoints);
		}
		if (nDim > 2)
		{
			zei = std::vector<double>(nIntPoints);
		}

		if (nGauss == 14)
		{
			for (int i = 0; i < nGauss; i++)
			{
				xii[i] = X14[i];
				eti[i] = Y14[i];
				zei[i] = Z14[i];
				wi[i] = (i < 8) ? Wa : Wb;
			}
		}
		else
		{
			int ip = 0;
			int n = nGauss - 1;
			switch (nDim)
			{
			case 1:
				for (int i = 0; i < nGauss; i++)
				{
					xii[ip] = X[n][i];
					wi[ip++] = W[n][i];
				}
				break;

			case 2:
				for (int i = 0; i < nGauss; i++)
				{
					for (int j = 0; j < nGauss; j++)
					{
						xii[ip] = X[n][i];
						eti[ip] = X[n][j];
						wi[ip++] = W[n][i] * W[n][j];
					}
				}
				break;

			case 3:
				for (int i = 0; i < nGauss; i++)
				{
					for (int j = 0; j < nGauss; j++)
					{
						for (int k = 0; k < nGauss; k++)
						{
							xii[ip] = X[n][i];
							eti[ip] = X[n][j];
							zei[ip] = X[n][k];
							wi[ip++] = W[n][i] * W[n][j] * W[n][k];
						}
					}
				}
				break;
			}
		}
	}
}
