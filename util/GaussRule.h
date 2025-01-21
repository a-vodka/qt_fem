#pragma once

#include <vector>
#include "rectangularvectors.h"

namespace util
{

	// Gauss integration rule
	class GaussRule
	{

		// Abscissas of the Gauss rule
	public:
		std::vector<double> xii, eti, zei;
		// Integration weights
		std::vector<double> wi;
		// Total namber of integration poins
		int nIntPoints = 0;

		// Abscissas and weights for 1, 2 and 3-point rules
	private:
		static const std::vector<std::vector<double>> X;
		static const std::vector<std::vector<double>> W;
		// Abscissas and weights for 14-point rule (3D)
		static constexpr double a = 0.7587869106393281;
			static constexpr double b = 0.7958224257542215;
		static const std::vector<double> X14;
		static const std::vector<double> Y14;
		static const std::vector<double> Z14;
		static constexpr double Wa = 0.3351800554016621;
			static constexpr double Wb = 0.8864265927977839;

		// Construct Gauss integration rule.
		// nGauss - number of Gauss points in each direction
		// (excluding 14-point rule),
		// nDim - number of dimensions
	public:
		GaussRule(int nGauss, int nDim);

	};

}
