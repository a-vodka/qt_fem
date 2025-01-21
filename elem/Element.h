#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "model/FeLoadData.h"
#include "rectangularvectors.h"

//JAVA TO C++ CONVERTER NOTE: Forward class declarations:
namespace model { class FeModel; }
namespace model { class FeLoad; }
namespace material { class Material; }
namespace elem { class StressContainer; }
namespace model { class ElemFaceLoad; }

namespace elem
{

	using namespace model;
	using namespace material;

	// Finite element
	class Element
	{

		// Finite element model
	public:
		static FeModel *fem;
		// Finite element load
        static FeLoadData *load;
		// Material of current element
		static Material *mat;
		// Element stiffness matrix
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: public static double kmat[][] = new double[60][60];
		static std::vector<std::vector<double>> kmat;
		// Element vector
		static std::vector<double> evec;
		// Element nodal coordinates
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: static double xy[][] = new double[20][3];
		static std::vector<std::vector<double>> xy;
		// Element nodal temperatures
		static std::vector<double> dtn;
		// Strain vector
		static std::vector<double> dstrain;

		// Element name
		std::wstring name;
		// Element material name
		std::wstring matName;
		// Element connectivities
		std::vector<int> ind;
		// Stress-strain storage
		std::vector<StressContainer*> str;

		// Implemented element types
	public:
		enum class elements
		{
//JAVA TO C++ CONVERTER TODO TASK: Enum value-specific class bodies are not converted by Java to C++ Converter:
//			quad8 {Element create() {return new ElementQuad2D();}},
//JAVA TO C++ CONVERTER TODO TASK: Enum value-specific class bodies are not converted by Java to C++ Converter:
//			hex20 {Element create() {return new ElementQuad3D();}};

//JAVA TO C++ CONVERTER TODO TASK: Enum methods are not converted by Java to C++ Converter:
//			abstract Element create();
		};

		// Construct new element
		// name - element name
	public:
		static Element *newElement(const std::wstring &name);

		// Constructor for an element.
		// name - element name;
		// nind - number of nodes;
		// nstress - number of stress points
		Element(const std::wstring &name, int nind, int nstress);

		// Compute element stiffness matrix kmat[][]
		virtual void stiffnessMatrix();

		// Compute element thermal vector (evec[])
		virtual void thermalVector();

		// Element nodal equivalent of distributed face load
		// (evec[])
		virtual int equivFaceLoad(ElemFaceLoad *surLd);

		// Nodal vector equivalent to stresses (evec[])
		virtual void equivStressVector();

		// Get local node numbers for element faces
		// returns elementFaces[nFaces][nNodesOnFace]
		virtual std::vector<std::vector<int>> getElemFaces();

		// Get strains at integration point (stress)
		// intPoint - integration point number (stress);
		// returns  strain vector [2*ndim]
		virtual std::vector<double> getStrainsAtIntPoint(int intPoint);

		// Get temperature at integration point (stress)
		// intPoint - integration point number (stress);
		// returns  temperature
		virtual double getTemperatureAtIntPoint(int intPoint);

		// Extrapolate quantity from integration points to nodes
		// fip [nInt][2*nDim] - values at integration points;
		// fn [nind][2*nDim] - values at nodes (out)
		virtual void extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn);

		// Set element connectivities
		// indel - connectivity numbers
		// nind - number of element nodes
		virtual void setElemConnectivities(std::vector<int> &indel, int nind);

		// Set element connectivities
		// indel - connectivity numbers
		virtual void setElemConnectivities(std::vector<int> &indel);

		// Set element material name
		// mat - material name
		virtual void setElemMaterial(const std::wstring &mat);

		// Set element nodal coordinates xy[nind][nDim]
		virtual void setElemXy();

		// Set nodal coordinates xy[nind][nDim] and
		//     temperatures dtn[nind]
		virtual void setElemXyT();

		// Assemble element vector.
		// elVector - element vector;
		// glVector - global vector (in/out)
		virtual void assembleElemVector(std::vector<double> &elVector, std::vector<double> &glVector);

		// Disassemble element vector (result in evec[]).
		// glVector - global vector
		virtual void disAssembleElemVector(std::vector<double> &glVector);

		// Returns element connectivities
		virtual std::vector<int> getElemConnectivities();

		//  Accumulate stresses and equivalent plastic strain
		virtual void accumulateStress();

	};

}
