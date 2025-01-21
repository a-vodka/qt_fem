#include "FE.h"

namespace fea
{
    int FE::main = 0;
    double FE::bigValue = 1.0e64;
    bool FE::tunedSolver = true;
    int FE::maxRow2D = 21;
    int FE::maxRow3D = 117;
    int FE::maxIterPcg = 10000;
    bool FE::epIntegrationTANGENT = false;
}
