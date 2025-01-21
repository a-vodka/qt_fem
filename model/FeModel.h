#pragma once

#include "FeModelData.h"
#include <string>
#include <vector>
#include <stdexcept>
#include "stringhelper.h"
#include "rectangularvectors.h"
#include "../solver/Solver.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace util { class FeScanner; }

namespace model
{

	using namespace elem;
	using namespace material;
	using namespace util;
	using namespace solver;


	// Description of the finite element model
	class FeModel : public FeModelData
	{

	private:
		static std::vector<int> elCon;
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: private static double box[][] = new double[2][3];
		static std::vector<std::vector<double>> box;
    public:
        std::list<Dof>::iterator *it;

		// Construct finite element model.
		// RD - data scanner, PR - print writer.
		virtual ~FeModel()
		{
			delete it;
		}
/*
		FeModel(FeScanner *RD, PrintWriter *PR);

		// Read data for a finite element model
		virtual void readData();

	private:
		void readDataFile(FeScanner *es);

		// Read element type, material and connectivities
		// for all elements
		void readElemData(FeScanner *es);

		// Read data for specified constrained displacements
		void readConstrDisplacements(FeScanner *es);

		// Create data for constrained displacements
		// specified inside a box
		void createBoxConstrDisplacements(FeScanner *es);
*/
	};

}
