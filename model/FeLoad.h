#pragma once

#include "FeLoadData.h"
#include <string>
#include <vector>
#include <list>
#include <stdexcept>
#include "stringhelper.h"
#include "rectangularvectors.h"
#include "../fea/FE.h"
#include "../elem/Element.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace model { class FeModel; }
namespace util { class FeScanner; }

namespace model
{

	using namespace fea;
	using namespace elem;
	using namespace util;


	// Load increment for the finite element model
	class FeLoad : public FeLoadData
	{

		// Finite element model
	private:
		static FeModel *fem;
    public:
        std::list<Dof>::iterator *itnf, *itsf;

		// Construct finite element load.
		// fem - finite element model
		virtual ~FeLoad()
		{
			delete itnf;
			delete itsf;
		}

		FeLoad(FeModel *fem);

		// Read data describing load increment.
		// returns  true if load data has been read
		virtual bool readData();

		// Read data fragment for load increment.
		// newLoad = true - beginning of new load,
		//         = false - continuation of load.
		// returns  true if load data has been read
	private:
		bool readDataFile(FeScanner *es, bool newLoad);


		// Read data for specified nodal forces
		void readNodalForces(FeScanner *es);

		// Read data for surface forces (element face loading):
		// direction, iel, nFaceNodes, faceNodes, forcesAtNodes.
		void readSurForces(FeScanner *es);

		// Create data for distributed surface load
		// specified inside a box
		void createBoxSurForces(FeScanner *es);

		// Assemble right-hand side of the global equation system
	public:
		virtual void assembleRHS();

	};

}
