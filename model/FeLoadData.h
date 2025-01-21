#pragma once

#include <string>
#include <vector>
#include <list>
#include <stdexcept>
#include <utility>
#include "rectangularvectors.h"
#include "../model/Dof.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace util { class FeScanner; }

namespace model
{

	using FeScanner = util::FeScanner;

	// Load data
	class FeLoadData
	{

	public:
		FeScanner *RD;
		static std::wstring loadStepName;
		// Load scale multiplier
		double scaleLoad = 0;
		// Relative residual norm tolerance
		static double residTolerance;
		// Maximum number of iterations (elastic-plastic problem)
		static int maxIterNumber;
        // Degrees of freedom with node forces
        std::list<Dof> nodForces;
        // Element face surface loads
        std::list<Dof> surForces;
		// Temperature increment
		static std::vector<double> dtemp;

		// Increment of force load
		static std::vector<double> dpLoad;
		// Total force load
		static std::vector<double> spLoad;
		// Increment of fictitious thermal loading
		static std::vector<double> dhLoad;
		// Displacement increment
		static std::vector<double> dDispl;
		// Total displacements
		static std::vector<double> sDispl;
		// Right-hand side of global equation system
		static std::vector<double> RHS;

		// Working arrays
		static std::vector<int> iw;
		static std::vector<double> dw;
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: static double[][] box = new double[2][3];
		static std::vector<std::vector<double>> box;

	public:
		enum class vars
		{
			loadstep,
			scaleload,
			residtolerance,
			maxiternumber,
			nodforce,
			surforce,
			boxsurforce,
			nodtemp,
			includefile,
			end
		};

		class varsHelper
		{
		private:
			static std::vector<std::pair<vars, std::wstring>> pairs()
			{
				return
				{
					{vars::loadstep, L"loadstep"},
					{vars::scaleload, L"scaleload"},
					{vars::residtolerance, L"residtolerance"},
					{vars::maxiternumber, L"maxiternumber"},
					{vars::nodforce, L"nodforce"},
					{vars::surforce, L"surforce"},
					{vars::boxsurforce, L"boxsurforce"},
					{vars::nodtemp, L"nodtemp"},
					{vars::includefile, L"includefile"},
					{vars::end, L"end"}
				};
			}

		public:
			static std::vector<vars> values()
			{
				std::vector<vars> temp;
				for (auto pair : pairs())
				{
					temp.push_back(pair.first);
				}
				return temp;
			}

			static std::wstring enumName(vars value)
			{
				for (auto pair : pairs())
				{
					if (pair.first == value)
						return pair.second;
				}

				throw std::runtime_error("Enum not found.");
			}

			static int ordinal(vars value)
			{
				std::vector<std::pair<vars, std::wstring>> temp = pairs();
				for (std::size_t i = 0; i < temp.size(); i++)
				{
					if (temp[i].first == value)
						return i;
				}

				throw std::runtime_error("Enum not found.");
			}

			static vars enumFromString(std::wstring value)
			{
				for (auto pair : pairs())
				{
					if (pair.second == value)
						return pair.first;
				}

				throw std::runtime_error("Enum not found.");
			}
		};


	public:
		virtual ~FeLoadData()
		{
            //delete RD;
		}

	};

}
