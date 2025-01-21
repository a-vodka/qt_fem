#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <list>
#include <stdexcept>
#include <utility>
#include "rectangularvectors.h"
#include "../material/Material.h"
#include "dof.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace util { class FeScanner; }
namespace elem { class Element; }

namespace model
{

	using namespace elem;
	using namespace util;


	// Finite element model data
	class FeModelData
	{

	public:
		static FeScanner *RD;
 //       static PrintWriter *PR;

		// Problem dimension =2/3
		int nDim = 3;
		// Number of degrees of freedom per node =2/3
		int nDf = 3;
		// Number of nodes
		int nNod = 0;
		// Number of elements
		int nEl = 0;
		// Number of degrees of freedom in the FE model
		int nEq = 0;
		// Elements
		std::vector<Element*> elems;
        // Materials
        std::unordered_map<std::wstring, material::Material*> materials;
		// Coordinates of nodes
	private:
		std::vector<double> xyz;
		// Constrained degrees of freedom
    public:
        std::list<Dof> defDs;
		bool thermalLoading = false;
		static std::wstring varName;

	public:
		enum class StrStates
		{
			plstrain,
			plstress,
			axisym,
			threed
		};

		class StrStatesHelper
		{
		private:
			static std::vector<std::pair<StrStates, std::wstring>> pairs()
			{
				return
				{
					{StrStates::plstrain, L"plstrain"},
					{StrStates::plstress, L"plstress"},
					{StrStates::axisym, L"axisym"},
					{StrStates::threed, L"threed"}
				};
			}

		public:
			static std::vector<StrStates> values()
			{
				std::vector<StrStates> temp;
				for (auto pair : pairs())
				{
					temp.push_back(pair.first);
				}
				return temp;
			}

			static std::wstring enumName(StrStates value)
			{
				for (auto pair : pairs())
				{
					if (pair.first == value)
						return pair.second;
				}

				throw std::runtime_error("Enum not found.");
			}

			static int ordinal(StrStates value)
			{
				std::vector<std::pair<StrStates, std::wstring>> temp = pairs();
				for (std::size_t i = 0; i < temp.size(); i++)
				{
					if (temp[i].first == value)
						return i;
				}

				throw std::runtime_error("Enum not found.");
			}

			static StrStates enumFromString(std::wstring value)
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
		static StrStates stressState;

	public:
		enum class PhysLaws
		{
			elastic,
			elplastic
		};

		class PhysLawsHelper
		{
		private:
			static std::vector<std::pair<PhysLaws, std::wstring>> pairs()
			{
				return
				{
					{PhysLaws::elastic, L"elastic"},
					{PhysLaws::elplastic, L"elplastic"}
				};
			}

		public:
			static std::vector<PhysLaws> values()
			{
				std::vector<PhysLaws> temp;
				for (auto pair : pairs())
				{
					temp.push_back(pair.first);
				}
				return temp;
			}

			static std::wstring enumName(PhysLaws value)
			{
				for (auto pair : pairs())
				{
					if (pair.first == value)
						return pair.second;
				}

				throw std::runtime_error("Enum not found.");
			}

			static int ordinal(PhysLaws value)
			{
				std::vector<std::pair<PhysLaws, std::wstring>> temp = pairs();
				for (std::size_t i = 0; i < temp.size(); i++)
				{
					if (temp[i].first == value)
						return i;
				}

				throw std::runtime_error("Enum not found.");
			}

			static PhysLaws enumFromString(std::wstring value)
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
		PhysLaws physLaw = PhysLaws::elastic;

		// Input data names
	public:
		enum class vars
		{
			nel,
			nnod,
			ndim,
			stressstate,
			physlaw,
			solver,
			elcon,
			nodcoord,
			material,
			constrdispl,
			boxconstrdispl,
			thermalloading,
			includefile,
			user,
			end
		};

		class varsHelper
		{
		private:
			static std::vector<std::pair<vars, std::wstring>> pairs()
			{
				return
				{
					{vars::nel, L"nel"},
					{vars::nnod, L"nnod"},
					{vars::ndim, L"ndim"},
					{vars::stressstate, L"stressstate"},
					{vars::physlaw, L"physlaw"},
					{vars::solver, L"solver"},
					{vars::elcon, L"elcon"},
					{vars::nodcoord, L"nodcoord"},
					{vars::material, L"material"},
					{vars::constrdispl, L"constrdispl"},
					{vars::boxconstrdispl, L"boxconstrdispl"},
					{vars::thermalloading, L"thermalloading"},
					{vars::includefile, L"includefile"},
					{vars::user, L"user"},
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


		// Allocation of nodal coordinate array
	public:
		virtual void newCoordArray();

		// Set coordinates of node
		virtual void setNodeCoords(int node, std::vector<double> &xyzn);

		// Set ith coordinates of node
		virtual void setNodeCoord(int node, int i, double v);

		// Get coordinates of node
		virtual std::vector<double> getNodeCoords(int node);

		// Get ith coordinate of node
		virtual double getNodeCoord(int node, int i);

	};

}
