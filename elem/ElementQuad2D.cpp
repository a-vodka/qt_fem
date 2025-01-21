#include "ElementQuad2D.h"
#include "../util/GaussRule.h"
#include "../model/FeModel.h"
#include "../material/Material.h"
#include "../util/UTIL.h"
#include "ShapeQuad2D.h"
#include "../model/ElemFaceLoad.h"
#include "../elem/Element.h"
#include "elem/StressContainer.h"

namespace elem
{
	using namespace model;
	using namespace material;
	using namespace util;
std::vector<std::vector<int>> ElementQuad2D::faceInd =
{
	{0,1,2},
	{2,3,4},
	{4,5,6},
	{6,7,0}
};
std::vector<double> ElementQuad2D::an(8);
std::vector<std::vector<double>> ElementQuad2D::dnxy = RectangularVectors::RectangularDoubleVector(8, 2);
std::vector<std::vector<double>> ElementQuad2D::bmat = RectangularVectors::RectangularDoubleVector(4, 16);
std::vector<std::vector<double>> ElementQuad2D::emat = RectangularVectors::RectangularDoubleVector(4, 4);
std::vector<double> ElementQuad2D::ept(4);
double ElementQuad2D::r = 0;
GaussRule *ElementQuad2D::gk = new GaussRule(2,2);
GaussRule *ElementQuad2D::gh = new GaussRule(3,2);
GaussRule *ElementQuad2D::gf = new GaussRule(3,1);
GaussRule *ElementQuad2D::gs = new GaussRule(2,2);

	ElementQuad2D::ElementQuad2D() : Element(L"quad8", 8, 4)
	{
	}

	void ElementQuad2D::stiffnessMatrix()
	{

		// Zeros to stiffness matrix kmat
		for (int i = 0; i < 16; i++)
		{
			for (int j = 0; j < 16; j++)
			{
							kmat[i][j] = 0;
			}
		}

		// ld = length of strain/stress vector (3 or 4)
		int ld = (FeModel::stressState == FeModel::StrStates::axisym) ? 4 : 3;
		// Material mat
        mat = fem->materials[matName];
		if (mat == nullptr)
		{
			UTIL::errorMsg(L"Element material name: " + matName);
		}
		mat->elasticityMatrix(emat);

		// Gauss integration loop
		for (int ip = 0; ip < gk->nIntPoints; ip++)
		{
			// Set displacement differentiation matrix bmat
			double det = setBmatrix(gk->xii[ip], gk->eti[ip]);
			double dv = det * gk->wi[ip];
			if (FeModel::stressState == FeModel::StrStates::axisym)
			{
				dv *= 2 * M_PI * r;
			}
			// Upper symmetrical part of the stiffness matrix
			for (int i = 0; i < 16; i++)
			{
				for (int j = i; j < 16; j++)
				{
					double s = 0;
					for (int k = 0; k < ld; k++)
					{
						for (int l = 0; l < ld; l++)
						{
							s += bmat[l][i] * emat[l][k] * bmat[k][j];
						}
					}
					kmat[i][j] += s * dv;
				}
			}
		}
	}

	double ElementQuad2D::setBmatrix(double xi, double et)
	{

		// Derivatives of shape functions
		double det = ShapeQuad2D::deriv(xi, et, ind, xy, dnxy);
		if (det <= 0)
		{
			UTIL::errorMsg(L"Negative/zero 8N element area");
		}
		if (FeModel::stressState == FeModel::StrStates::axisym)
		{
			ShapeQuad2D::shape(xi, et, ind, an);
			r = 0;
			for (int i = 0; i < 8; i++)
			{
				r += an[i] * xy[i][0];
			}
		}
		// Eight blocks of the displacement differentiation
		//   matrix
		for (int ib = 0; ib < 8; ib++)
		{
			bmat[0][2 * ib] = dnxy[ib][0];
			bmat[0][2 * ib + 1] = 0.0;
			bmat[1][2 * ib] = 0.0;
			bmat[1][2 * ib + 1] = dnxy[ib][1];
			bmat[2][2 * ib] = dnxy[ib][1];
			bmat[2][2 * ib + 1] = dnxy[ib][0];
			if (FeModel::stressState == FeModel::StrStates::axisym)
			{
				bmat[3][2 * ib] = an[ib] / r;
				bmat[3][2 * ib + 1] = 0.0;
			}
		}
		return det;
	}

	void ElementQuad2D::thermalVector()
	{

		// Zeros to thermal vector evec
		for (int i = 0; i < 16; i++)
		{
			evec[i] = 0.0;
		}
		int ld = (FeModel::stressState == FeModel::StrStates::axisym) ? 4 : 3;
		// Material mat
		mat = static_cast<Material*>(fem->materials[matName]);
		mat->elasticityMatrix(emat);
		double alpha = mat->getAlpha();
		double nu = mat->getNu();

		// Gauss integration loop
		for (int ip = 0; ip < gh->nIntPoints; ip++)
		{
			// Set displacement differentiation matrix bmat
			double det = setBmatrix(gh->xii[ip], gh->eti[ip]);
			// Shape functions an
			ShapeQuad2D::shape(gh->xii[ip], gh->eti[ip], ind, an);
			double t = 0;
			for (int i = 0; i < 8; i++)
			{
				t += an[i] * dtn[i];
			}
			double dv = det * gh->wi[ip];
			if (FeModel::stressState == FeModel::StrStates::axisym)
			{
				dv *= 2 * M_PI * r;
			}
			ept[0] = alpha * t;
			if (FeModel::stressState == FeModel::StrStates::plstrain)
			{
					ept[0] *= (1 + nu);
			}
			ept[1] = ept[0];
			ept[2] = 0.0;
			ept[3] = ept[0];

			for (int i = 0; i < 16; i++)
			{
				double s = 0;
				for (int j = 0; j < ld; j++)
				{
					for (int k = 0; k < ld; k++)
					{
						s += bmat[k][i] * emat[j][k] * ept[j];
					}
				}
				evec[i] += s * dv;
			}
		}
	}

	int ElementQuad2D::equivFaceLoad(ElemFaceLoad *surLd)
	{
		// Shape functons
		std::vector<double> an(3);
		// Derivatives of shape functions
		std::vector<double> xin(3);

		for (int i = 0; i < 16; i++)
		{
			evec[i] = 0.0;
		}
		int loadedFace = surLd->rearrange(faceInd, ind);
		if (loadedFace == -1)
		{
			return -1;
		}

		// Gauss integration loop
		for (int ip = 0; ip < gf->nIntPoints; ip++)
		{
			ShapeQuad2D::shapeDerivFace(gf->xii[ip], surLd->faceNodes[1], an, xin);
			double p = r = 0;
			double xs = 0;
			double ys = 0;
			for (int i = 0; i < 3; i++)
			{
				p += an[i] * surLd->forceAtNodes[i];
				int j = faceInd[loadedFace][i];
				r += an[i] * xy[j][0];
				xs += xin[i] * xy[j][0];
				ys += xin[i] * xy[j][1];
			}
			double dl = std::sqrt(xs * xs + ys * ys);
			double ds = dl;
			if (FeModel::stressState == FeModel::StrStates::axisym)
			{
				 ds *= 2 * M_PI * r;
			}
			double p1, p2;
			// direction=0 - normal load, =1,2 - along axes x,y
			if (surLd->direction == 0 && ds > 0.0)
			{
				p1 = p * ys / dl;
				p2 = -p * xs / dl;
			}
			else if (surLd->direction == 1)
			{
				p1 = p;
				p2 = 0;
			}
			else
			{
				p1 = 0;
				p2 = p;
			}

			for (int i = 0; i < 3; i++)
			{
				int j = faceInd[loadedFace][i];
				evec[2 * j] += an[i] * p1 * ds * gf->wi[ip];
				evec[2 * j + 1] += an[i] * p2 * ds * gf->wi[ip];
			}
		}
		return loadedFace;
	}

	void ElementQuad2D::equivStressVector()
	{

		for (int i = 0; i < 16; i++)
		{
			evec[i] = 0.0;
		}
		int ld = (FeModel::stressState == FeModel::StrStates::axisym) ? 4 : 3;

		for (int ip = 0; ip < gs->nIntPoints; ip++)
		{
			// Accumulated stress
			std::vector<double> s(4);
			for (int i = 0; i < 4; i++)
			{
				s[i] = str[ip]->sStress[i] + str[ip]->dStress[i];
			}
			// Set displacement differentiation matrix bmat
			double det = setBmatrix(gs->xii[ip], gs->eti[ip]);
			double dv = det * gs->wi[ip];
			if (FeModel::stressState == FeModel::StrStates::axisym)
			{
				dv *= 2 * M_PI * r;
			}

			for (int i = 0; i < 16; i++)
			{
				double a = 0;
				for (int j = 0; j < ld; j++)
				{
					a += bmat[j][i] * s[j];
				}
				evec[i] -= a * dv;
			}
		}
	}

	void ElementQuad2D::extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn)
	{
		constexpr double A = 1 + 0.5 * std::sqrt(3.0), B = -0.5, C = 1 - 0.5 * std::sqrt(3.0);
		// Extrapolation matrix
		const std::vector<std::vector<double>> lim =
		{
			std::vector<double> {A, B, B, C},
			std::vector<double> {B, C, A, B},
			std::vector<double> {C, B, B, A},
			std::vector<double> {B, A, C, B}
		};

		for (int i = 1; i < 8; i += 2)
		{
			for (int j = 0; j < 4; j++)
			{
				fn[i][j] = 0;
			}
		}

		for (int corner = 0; corner < 4; corner++)
		{
			int n = (corner == 0) ? 7 : 2 * corner - 1;
			for (int k = 0; k < 4; k++)
			{
				double c = 0.0;
				for (int ip = 0 ; ip < 4; ip++)
				{
					c += lim[corner][ip] * fip[ip][k];
				}
				fn[2 * corner][k] = c; // corner node
				fn[n][k] += 0.5 * c;
				fn[2 * corner + 1][k] += 0.5 * c;
			}
		}
	}

	std::vector<std::vector<int>> ElementQuad2D::getElemFaces()
	{
		return faceInd;
	}

	std::vector<double> ElementQuad2D::getStrainsAtIntPoint(int ip)
	{

		// Set displacement differentiation matrix bmat
		setBmatrix(gs->xii[ip], gs->eti[ip]);
		std::vector<double> strain(4);
		for (int i = 0; i < 4; i++)
		{
			strain[i] = 0;
			for (int j = 0; j < 16; j++)
			{
				strain[i] += bmat[i][j] * evec[j];
			}
		}
		return strain;
	}

	double ElementQuad2D::getTemperatureAtIntPoint(int ip)
	{
		ShapeQuad2D::shape(gs->xii[ip], gs->eti[ip], ind, an);
		double t = 0;
		for (int i = 0; i < 8; i++)
		{
			t += an[i] * dtn[i];
		}
		return t;
	}
}
