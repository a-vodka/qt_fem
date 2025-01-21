#include "ElemFaceLoad.h"

namespace model
{

	ElemFaceLoad::ElemFaceLoad(int iel, int nFaceNodes, int direction, std::vector<int> &faceNodes, std::vector<double> &forceAtNodes)
	{
		this->iel = iel;
		this->direction = direction;
		this->faceNodes = std::vector<int>(nFaceNodes);
		this->forceAtNodes = std::vector<double>(nFaceNodes);

		for (int i = 0; i < nFaceNodes; i++)
		{
			this->faceNodes[i] = faceNodes[i];
			this->forceAtNodes[i] = forceAtNodes[i];
		}
	}

	ElemFaceLoad::ElemFaceLoad(int iel, int nFaceNodes, int direction, std::vector<int> &faceNodes, double force)
	{
		this->iel = iel;
		this->direction = direction;
		this->faceNodes = std::vector<int>(nFaceNodes);
		this->forceAtNodes = std::vector<double>(nFaceNodes);

		for (int i = 0; i < nFaceNodes; i++)
		{
			this->faceNodes[i] = faceNodes[i];
			this->forceAtNodes[i] = force;
		}
	}

	int ElemFaceLoad::rearrange(std::vector<std::vector<int>> &faces, std::vector<int> &ind)
	{

		std::vector<int> perm(8);
		std::vector<double> fw(8);
		int loadedFace = -1;

        for (std::size_t iface = 0; iface < faces.size(); iface++)
        {
            int nNodes = faces[iface].size();
			for (int inod = 0; inod < nNodes; inod++)
			{
				perm[inod] = -1;
			}
			for (int inod = 0; inod < nNodes; inod++)
			{
				int iGlob = ind[faces[iface][inod]];
				if (iGlob > 0)
				{
					bool EQ = false;
					int i;
					for (i = 0; i < nNodes; i++)
					{
						if (faceNodes[i] == iGlob)
						{
							EQ = true;
							break;
						}
					}
					if (!EQ)
					{
						goto FACEContinue;
					}
					perm[inod] = i;
				}
			}
			loadedFace = iface;
			for (int inod = 0; inod < nNodes; inod++)
			{
				faceNodes[inod] = ind[faces[iface][inod]];
				fw[inod] = forceAtNodes[inod];
			}
			for (int inod = 0; inod < nNodes; inod++)
			{
				forceAtNodes[inod] = (perm[inod] == -1) ? 0.0 : fw[perm[inod]];
			}
			FACEContinue:;
		}
		FACEBreak:
		return loadedFace;
	}
}
