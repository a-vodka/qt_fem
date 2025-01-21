#pragma once

#include <string>
#include <stdexcept>
#include "rectangularvectors.h"

namespace util
{


	//  Finite element printer to file
	class FePrintWriter
	{
	public:
        //PrintWriter *PR;
		virtual ~FePrintWriter()
		{
            //delete PR;
		}

        //virtual PrintWriter *getPrinter(const std::wstring &fileOut);

	};

}
