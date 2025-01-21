#include "FeLoadData.h"
//#include "../util/FeScanner.h"

namespace model
{
//	using FeScanner = util::FeScanner;
    std::wstring FeLoadData::loadStepName;
    double FeLoadData::residTolerance = 0.01;
    int FeLoadData::maxIterNumber = 100;
    std::vector<double> FeLoadData::dtemp;
    std::vector<double> FeLoadData::dpLoad;
    std::vector<double> FeLoadData::spLoad;
    std::vector<double> FeLoadData::dhLoad;
    std::vector<double> FeLoadData::dDispl;
    std::vector<double> FeLoadData::sDispl;
    std::vector<double> FeLoadData::RHS;
    std::vector<int> FeLoadData::iw(8);
    std::vector<double> FeLoadData::dw(8);
    std::vector<std::vector<double>> FeLoadData::box = RectangularVectors::RectangularDoubleVector(2, 3);
}
