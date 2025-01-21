#include "ElementQuad3D.h"
#include "../util/GaussRule.h"
#include "../material/Material.h"
#include "../util/UTIL.h"
#include "ShapeQuad3D.h"
#include "../model/ElemFaceLoad.h"
#include "../elem/Element.h"
#include "elem/StressContainer.h"
#include "model/FeModel.h"

namespace elem
{
	using namespace model;
	using namespace material;
	using namespace util;
std::vector<std::vector<int>> ElementQuad3D::faceInd =
{
	{0, 8,12,19,18,11, 6, 7},
	{2, 3, 4,10,16,15,14, 9},
	{0, 1, 2, 9,14,13,12, 8},
	{4, 5, 6,11,18,17,16,10},
	{0, 7, 6, 5, 4, 3, 2, 1},
	{12,13,14,15,16,17,18,19}
};
std::vector<double> ElementQuad3D::an(20);
std::vector<std::vector<double>> ElementQuad3D::dnxy = RectangularVectors::RectangularDoubleVector(20, 3);
GaussRule *ElementQuad3D::gk = new GaussRule(14,3);
GaussRule *ElementQuad3D::gh = new GaussRule(3,3);
GaussRule *ElementQuad3D::gf = new GaussRule(3,2);
GaussRule *ElementQuad3D::gs = new GaussRule(2,3);

	ElementQuad3D::ElementQuad3D() : Element(L"hex20", 20, 8)
	{
	}

	void ElementQuad3D::stiffnessMatrix()
	{

		for (int i = 0; i < 60; i++)
		{
			for (int j = i; j < 60; j++)
			{
				kmat[i][j] = 0.0;
			}
		}
		// Material mat
		mat = static_cast<Material*>(fem->materials[matName]);
		if (mat == nullptr)
		{
			UTIL::errorMsg(L"Element material name: " + matName);
		}
		double lambda = mat->getLambda();
		double mu = mat->getMu();
		double beta = lambda + 2 * mu;

		for (int ip = 0; ip < gk->nIntPoints; ip++)
		{
			double det = ShapeQuad3D::deriv(gk->xii[ip], gk->eti[ip], gk->zei[ip], ind, xy, dnxy);
			double dv = det * gk->wi[ip];
			// Upper symmetrical part of the matrix by rows
			for (int i = 0; i < 20; i++)
			{ // i = row
				// dNi/dx, dNi/dy, dNi/dz
				double dix = dnxy[i][0];
				double diy = dnxy[i][1];
				double diz = dnxy[i][2];
				for (int j = i; j < 20; j++)
				{ // j = column
					// dNj/dx, dNj/dy, dNj/dz
					double djx = dnxy[j][0];
					double djy = dnxy[j][1];
					double djz = dnxy[j][2];

					kmat[i * 3][j * 3] += (beta * dix * djx + mu * (diy * djy + diz * djz)) * dv;
					kmat[i * 3][j * 3 + 1] += (lambda * dix * djy + mu * diy * djx) * dv;
					kmat[i * 3][j * 3 + 2] += (lambda * dix * djz + mu * diz * djx) * dv;

					if (j > i)
					{
						kmat[i * 3 + 1][j * 3] += (lambda * diy * djx + mu * dix * djy) * dv;
					}
					kmat[i * 3 + 1][j * 3 + 1] += (beta * diy * djy + mu * (diz * djz + dix * djx)) * dv;
					kmat[i * 3 + 1][j * 3 + 2] += (lambda * diy * djz + mu * diz * djy) * dv;

					if (j > i)
					{
						kmat[i * 3 + 2][j * 3] += (lambda * diz * djx + mu * dix * djz) * dv;
						kmat[i * 3 + 2][j * 3 + 1] += (lambda * diz * djy + mu * diy * djz) * dv;
					}
					kmat[i * 3 + 2][j * 3 + 2] += (beta * diz * djz + mu * (dix * djx + diy * djy)) * dv;
				}
			}
		}
	}

	void ElementQuad3D::thermalVector()
	{

		for (int i = 0; i < 60; i++)
		{
			evec[i] = 0.0;
		}

		mat = static_cast<Material*>(fem->materials[matName]);
		double alpha = mat->getAlpha();
		double lambda = mat->getLambda();
		double mu = mat->getMu();
		double g = 3 * lambda + 2 * mu;

		for (int ip = 0; ip < gh->nIntPoints; ip++)
		{
			ShapeQuad3D::shape(gh->xii[ip], gh->eti[ip], gh->zei[ip], ind, an);
			// Temperature at integration point
			double t = 0;
			for (int i = 0; i < 20; i++)
			{
				t += an[i] * dtn[i];
			}
			double det = ShapeQuad3D::deriv(gh->xii[ip], gh->eti[ip], gh->zei[ip], ind, xy, dnxy);
			double dv = g * alpha * t * det * gh->wi[ip];
			for (int i = 0; i < 20; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					evec[i * 3 + j] += dnxy[i][j] * dv;
				}
			}
		}
	}

	int ElementQuad3D::equivFaceLoad(ElemFaceLoad *surLd)
	{
		// Shape functons
		std::vector<double> an(8);
		// Derivatives of shape functions
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: double xin[][] = new double[8][2];
		std::vector<std::vector<double>> xin = RectangularVectors::RectangularDoubleVector(8, 2);
		// Tangent vectors along xi and eta
//JAVA TO C++ CONVERTER NOTE: The following call to the 'RectangularVectors' helper class reproduces the rectangular array initialization that is automatic in Java:
//ORIGINAL LINE: double e[][] = new double[2][3];
		std::vector<std::vector<double>> e = RectangularVectors::RectangularDoubleVector(2, 3);
		// Normal vector
		std::vector<double> g(3);
		std::vector<double> ps(3);

		for (int i = 0; i < 60; i++)
		{
			evec[i] = 0.0;
		}

		int loadedFace = surLd->rearrange(faceInd, ind);
		if (loadedFace == -1)
		{
			return -1;
		}

		for (int ip = 0; ip < gf->nIntPoints; ip++)
		{
			ShapeQuad3D::shapeDerivFace(gf->xii[ip], gf->eti[ip], ind, an, xin);
			double p = 0.0;
			for (int i = 0; i < 8; i++)
			{
				p += an[i] * surLd->forceAtNodes[i];
			}
			// Tangent vectors
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					double s = 0;
					for (int k = 0; k < 8; k++)
					{
						s += xin[k][i] * xy[faceInd[loadedFace][k]][j];
					}
					e[i][j] = s;
				}
			}
			 // Normal vector g
			 g[0] = (e[0][1] * e[1][2] - e[1][1] * e[0][2]);
			 g[1] = (e[0][2] * e[1][0] - e[1][2] * e[0][0]);
			 g[2] = (e[0][0] * e[1][1] - e[1][0] * e[0][1]);

			// Element of surface ds
			double ds = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
			if (ds <= 0)
			{
				UTIL::errorMsg(L"Negative/zero element face");
			}
			// Surface load components ps:
			// direction=0 - normal, x=1, y=2, z=3
			if (surLd->direction == 0)
			{
				for (int i = 0; i < 3; i++)
				{
					ps[i] = p * g[i] / ds;
				}
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					ps[i] = 0;
				}
				ps[surLd->direction - 1] = p;
			}
			for (int i = 0; i < 8; i++)
			{
				int k = faceInd[loadedFace][i];
				for (int j = 0; j < 3; j++)
				{
					evec[3 * k + j] += an[i] * ps[j] * ds * gf->wi[ip];
				}
			}
		}
		return loadedFace;
	}

	void ElementQuad3D::equivStressVector()
	{

		for (int i = 0; i < 60; i++)
		{
			evec[i] = 0.0;
		}

		for (int ip = 0; ip < gs->nIntPoints; ip++)
		{
			// Accumulated stress  s
			std::vector<double> s(6);
			for (int i = 0; i < 6; i++)
			{
				s[i] = str[ip]->sStress[i] + str[ip]->dStress[i];
			}
			double det = ShapeQuad3D::deriv(gs->xii[ip], gs->eti[ip], gs->zei[ip], ind, xy, dnxy);
			double dv = det * gs->wi[ip];

			for (int i = 0; i < 20; i++)
			{
				double a0 = dnxy[i][0];
				double a1 = dnxy[i][1];
				double a2 = dnxy[i][2];
				evec[i * 3] -= (a0 * s[0] + a1 * s[3] + a2 * s[5]) * dv;
				evec[i * 3 + 1] -= (a1 * s[1] + a0 * s[3] + a2 * s[4]) * dv;
				evec[i * 3 + 2] -= (a2 * s[2] + a1 * s[4] + a0 * s[5]) * dv;
			}
		}
	}

	void ElementQuad3D::extrapolateToNodes(std::vector<std::vector<double>> &fip, std::vector<std::vector<double>> &fn)
	{
		// Vertices
		const std::vector<int> vn = {0, 2, 4, 6, 12, 14, 16, 18};
		 // Midside nodes
		const std::vector<int> mn = {8, 9, 10, 11, 8, 9, 10, 11};
		// Extrapolation matrix
		constexpr double A = 0.25 * (5 + 3 * std::sqrt(3.0)), B = -0.25 * (std::sqrt(3.0) + 1), C = 0.25 * (std::sqrt(3.0) - 1), D = 0.25 * (5 - 3 * std::sqrt(3.0));
		const std::vector<std::vector<double>> lim =
		{
			std::vector<double> {A, B, B, C, B, C, C, D},
			std::vector<double> {B, C, C, D, A, B, B, C},
			std::vector<double> {C, D, B, C, B, C, A, B},
			std::vector<double> {B, C, A, B, C, D, B, C},
			std::vector<double> {B, A, C, B, C, B, D, C},
			std::vector<double> {C, B, D, C, B, A, C, B},
			std::vector<double> {D, C, C, B, C, B, B, A},
			std::vector<double> {C, B, B, A, D, C, C, B}
		};

		for (int i = 0; i < 20; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				fn[i][j] = 0;
			}
		}

		for (int vertex = 0; vertex < 8; vertex++)
		{
			int i = vn[vertex]; // node at vertex
			int im = i - 1;
			if (i == 0)
			{
				im = 7;
			}
			if (i == 12)
			{
				im = 19;
			}
			for (int k = 0; k < 6; k++)
			{
				double c = 0.0;
				for (int j = 0 ; j < 8; j++)
				{
					c += fip[j][k] * lim[vertex][j];
				}
				fn[i][k] = c;
				fn[im][k] += 0.5 * c;
				fn[i + 1][k] += 0.5 * c;
				fn[mn[vertex]][k] += 0.5 * c;
			}
		}
	}

	std::vector<std::vector<int>> ElementQuad3D::getElemFaces()
	{
		return faceInd;
	}

	std::vector<double> ElementQuad3D::getStrainsAtIntPoint(int ip)
	{

		// Derivatives of shape functions
		ShapeQuad3D::deriv(gs->xii[ip], gs->eti[ip], gs->zei[ip], ind, xy, dnxy);

		// Derivatives of displacements
		double dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz;
		dux = duy = duz = dvx = dvy = dvz = dwx = dwy = dwz = 0;
		for (int i = 0; i < 20; i++)
		{
			double dnx = dnxy[i][0];
			double dny = dnxy[i][1];
			double dnz = dnxy[i][2];
			double u = evec[3 * i];
			double v = evec[3 * i + 1];
			double w = evec[3 * i + 2];
			dux += dnx * u;
			duy += dny * u;
			duz += dnz * u;
			dvx += dnx * v;
			dvy += dny * v;
			dvz += dnz * v;
			dwx += dnx * w;
			dwy += dny * w;
			dwz += dnz * w;
		}
		// Strains
		std::vector<double> strain(6);
		strain[0] = dux;
		strain[1] = dvy;
		strain[2] = dwz;
		strain[3] = duy + dvx;
		strain[4] = dvz + dwy;
		strain[5] = duz + dwx;
		return strain;
	}

	double ElementQuad3D::getTemperatureAtIntPoint(int ip)
	{

		ShapeQuad3D::shape(gs->xii[ip], gs->eti[ip], gs->zei[ip], ind, an);
		double t = 0;
		for (int i = 0; i < 20; i++)
		{
			t += an[i] * dtn[i];
		}
		return t;
	}
}
