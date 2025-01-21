#pragma once

#include <string>
#include <iostream>
#include "rectangularvectors.h"

namespace util
{
    class PrintWriter {

    };

	// Miscellaneous static methods
	class UTIL
	{

		// Print date and time.
		// PR - PrintWriter for listing file
	public:
		static void printDate(PrintWriter *PR);

		// Print error message and exit.
		// message - error message that is printed.
		static void errorMsg(const std::wstring &message);

		// Transform text direction into integer.
		// s - direction x/y/z/n.
		// returns  integer direction 1/2/3/0, error: -1.
		static int direction(const std::wstring &s);

	};

}
