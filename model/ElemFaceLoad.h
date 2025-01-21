#pragma once

#include <vector>
#include "rectangularvectors.h"

namespace model
{

	// Element face load
	class ElemFaceLoad
	{
		// Element number (start with 0)
	public:
		int iel = 0;
		// Direction: 1-x, 2-y, 3-z, 0-normal
		int direction = 0;
		std::vector<int> faceNodes;
		std::vector<double> forceAtNodes;

		ElemFaceLoad(int iel, int nFaceNodes, int direction, std::vector<int> &faceNodes, std::vector<double> &forceAtNodes);

		ElemFaceLoad(int iel, int nFaceNodes, int direction, std::vector<int> &faceNodes, double force);

		// Rearrange surface load (faceNodes[] and ForcesAtNodes[])
		// according to order in element faces.
		// faces - local numbers (from zero) of element faces,
		// ind - element connectivities.
		// returns  loaded face number or -1 if no match
		//  between ind[] and load data.
		virtual int rearrange(std::vector<std::vector<int>> &faces, std::vector<int> &ind);

	};
}
