#include <iostream>

#include "model/FeLoad.h"
#include "model/FeLoadData.h"
#include "model/FeModel.h"
#include "model/FeStress.h"

using namespace std;

using namespace elem;
using namespace model;
using namespace solver;
using namespace util;

int main()
{

    FeModel *fem = new FeModel();
    Element::fem = fem;

    Material* mat = Material::newMaterial(L"elastic", L"threed");
    double e = 2.1e11;
    double nu = 0.3;
    double alpha = 0;
    mat->setElasticProp(e, nu, alpha);
    fem->materials[L"1"]=mat;

    const int nNod = 20;
    const int nDim = 3;


    fem->nNod = nNod;
    fem->nDim = nDim;
    fem->nEl = 1;
    fem->nEq = fem->nNod * fem->nDim;
    fem->newCoordArray();

    double node[nNod][nDim] =
    {
        {0.0, 0.0 ,0.0}, //1
        {0.5, 0.0, 0.0}, //2
        {1.0, 0.0, 0.0}, //3
        {1.0, 0.5, 0.0}, //4
        {1.0, 1.0, 0.0}, //5
        {0.5, 1.0 ,0.0}, //6
        {0.0, 1.0 ,0.0}, //7
        {0.0, 0.5 ,0.0}, //8
        {0.0, 0.0 ,0.5}, //9
        {1.0, 0.0 ,0.5}, //10
        {1.0, 1.0 ,0.5}, //11
        {0.0, 1.0 ,0.5}, //12
        {0.0, 0.0 ,1.0}, //13
        {0.5, 0.0, 1.0}, //14
        {1.0, 0.0, 1.0}, //15
        {1.0, 0.5, 1.0}, //16
        {1.0, 1.0, 1.0}, //17
        {0.5, 1.0 ,1.0}, //18
        {0.0, 1.0 ,1.0}, //19
        {0.0, 0.5 ,1.0}  //20
    };

    // Nodal coordinates
    for (int i = 0; i < nNod; i++)
    {
        for (int j = 0; j < nDim; j++)
        {
            fem->setNodeCoord(i, j, node[i][j]);
        }
    }

    fem->elems.resize(fem->nEl);
    for (int iel = 0; iel < fem->nEl; iel++)
    {
        // Element type
        std::wstring s = L"hex20";
        fem->elems[iel] = Element::newElement(s);
        // Element material
        std::wstring elMat = L"1";
        fem->elems[iel]->setElemMaterial(elMat);
        // Element connectivities
        int nind = fem->elems[iel]->ind.size();
        vector<int> elCon;
        for (int l = 1; l<= 20; l++)
            elCon.push_back(l);
        fem->elems[iel]->setElemConnectivities(elCon,nind);
    }


    Solver *solver = Solver::newSolver(fem);
    solver->assembleGSM();

    wcout<<L"Memory for global matrix: " << Solver::lengthOfGSM * 8.0e-6<<" MB\n" << Solver::lengthOfGSM <<endl;



    FeLoadData *load = new FeLoadData();
    Element::load = load;



    load->nodForces.push_back(Dof(fem->nDim*(3-1)+UTIL::direction(L"x"), 1e6));
    load->nodForces.push_back(Dof(fem->nDim*(5-1)+UTIL::direction(L"x"), 1e6));
    load->nodForces.push_back(Dof(fem->nDim*(15-1)+UTIL::direction(L"x"), 1e6));
    load->nodForces.push_back(Dof(fem->nDim*(17-1)+UTIL::direction(L"x"), 1e6));


    fem->defDs.push_back(Dof(fem->nDim*(1-1)+UTIL::direction(L"x"), 0.0));
    fem->defDs.push_back(Dof(fem->nDim*(7-1)+UTIL::direction(L"x"), 0.0));
    fem->defDs.push_back(Dof(fem->nDim*(13-1)+UTIL::direction(L"x"), 0.0));
    fem->defDs.push_back(Dof(fem->nDim*(19-1)+UTIL::direction(L"x"), 0.0));
    fem->defDs.push_back(Dof(fem->nDim*(1-1)+UTIL::direction(L"y"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(7-1)+UTIL::direction(L"y"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(13-1)+UTIL::direction(L"y"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(19-1)+UTIL::direction(L"y"), 0.0));
    fem->defDs.push_back(Dof(fem->nDim*(1-1)+UTIL::direction(L"z"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(7-1)+UTIL::direction(L"z"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(13-1)+UTIL::direction(L"z"), 0.0));
    //fem->defDs.push_back(Dof(fem->nDim*(19-1)+UTIL::direction(L"z"), 0.0));
    /*load->nodForces.push_back(Dof(fem->nDim*(1-1)+UTIL::direction(L"x"), -1e6));
    load->nodForces.push_back(Dof(fem->nDim*(7-1)+UTIL::direction(L"x"), -1e6));
    load->nodForces.push_back(Dof(fem->nDim*(13-1)+UTIL::direction(L"x"), -1e6));
    load->nodForces.push_back(Dof(fem->nDim*(19-1)+UTIL::direction(L"x"), -1e6));*/

    FeStress *stress = new FeStress(fem);
    load->dpLoad.resize(fem->nEq);
    load->RHS.resize(fem->nEq);
    load->dhLoad.resize(fem->nEq);
    load->dDispl.resize(fem->nEq);
    load->spLoad.resize(fem->nEq);
    load->sDispl.resize(fem->nEq);
    for (auto i = load->nodForces.begin(); i!=load->nodForces.end(); i++)
    {
//        wcout << i->dofNum << '\t' << i->value << endl << flush;
        load->dpLoad[i->dofNum - 1] = i->value;
    }

    // Right-hand side = actual load + fictitious load
    for (int i = 0; i < fem->nEq; i++)
    {
        load->RHS[i] = load->dpLoad[i] + load->dhLoad[i];
    }


    for (auto i = fem->defDs.begin();i!=fem->defDs.end();i++)
    {
        load->RHS[i->dofNum - 1] = FE::bigValue * i->value;
    }

    //while (load->readData())
    {
        //load->assembleRHS();
        int iter = 0;
        // Equilibrium iterations
        do
        {
            iter++;
            int its = solver->solve(FeLoad::RHS);
            if (its > 0)
            {
                wcout<<L"Solver: " << its << L"iterations\n";
            }
            stress->computeIncrement();
        } while (!stress->equilibrium(iter));

        stress->accumulate();
        stress->writeResults();
        wcout << L"Loadstep " << FeLoad::loadStepName;
        if (iter > 1)
        {
            wcout << iter <<L" %5d iterations, Relative residual norm = " << FeStress::relResidNorm;
        }
        wcout << endl;
    }

    //PR->printf(L"\nSolution time = %10.2f s\n", (System::currentTimeMillis() - t0) * 0.001);



    cout << "Hello World!" << endl;
    return 0;
}
