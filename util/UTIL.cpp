#include "UTIL.h"

namespace util
{

	void UTIL::printDate(PrintWriter *PR)
	{
        /*

		Calendar *c = new GregorianCalendar();

		PR->printf(L"Date: %d-%02d-%02d  Time: %02d:%02d:%02d\n", c->get(Calendar::YEAR), c->get(Calendar::MONTH) + 1, c->get(Calendar::DATE), c->get(Calendar::HOUR_OF_DAY), c->get(Calendar::MINUTE),c->get(Calendar::SECOND));

        delete c;*/
	}

	void UTIL::errorMsg(const std::wstring &message)
	{
			std::wcout << L"=== ERROR: " << message << std::endl;
			exit(1);
	}

	int UTIL::direction(const std::wstring &s)
	{
		if (s == L"x")
		{
			return 1;
		}
		else if (s == L"y")
		{
			return 2;
		}
		else if (s == L"z")
		{
			return 3;
		}
		else if (s == L"n")
		{
			return 0;
		}
		else
		{
			return -1;
		}
	}
}
