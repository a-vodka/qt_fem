#include "FeScanner.h"
#include "UTIL.h"
#include "../model/Dof.h"

namespace util
{
	using Dof = model::Dof;

	FeScanner::FeScanner(const std::wstring &fileIn)
	{

		try
		{
			File tempVar(fileIn);
			es = new Scanner(&tempVar);
		}
		catch (const std::runtime_error &e)
		{
			UTIL::errorMsg(L"Input file not found: " + fileIn);
		}
		es->useDelimiter(L"\\s*=\\s*|\\s+");

	}

	bool FeScanner::hasNext()
	{
		return es->hasNext();
	}

	bool FeScanner::hasNextDouble()
	{
		return es->hasNextDouble();
	}

	std::wstring FeScanner::next()
	{
		return es->next();
	}

	double FeScanner::nextDouble()
	{
		return es->nextDouble();
	}

	int FeScanner::readInt()
	{
		if (!es->hasNextInt())
		{
			UTIL::errorMsg(L"Expected integer. Instead: " + es->next());
		}
		return es->nextInt();
	}

	double FeScanner::readDouble()
	{
		if (!es->hasNextDouble())
		{
			UTIL::errorMsg(L"Expected double. Instead: " + es->next());
		}
		return es->nextDouble();
	}

	void FeScanner::nextLine()
	{
		es->nextLine();
	}

	void FeScanner::moveAfterLineWithWord(const std::wstring &word)
	{

		while (es->hasNext())
		{
			std::wstring varname = es->next()->toLowerCase();
			if (varname == L"#")
			{
				es->nextLine();
											  continue;
			}
			if (varname == word)
			{
				es->nextLine();
				return;
			}
			es++;
		}
		UTIL::errorMsg(L"moveAfterLineWithWord cannot find: " + word);
	}

	ListIterator *FeScanner::readNumberList(ListIterator *it, int dir, int ndim, double sValue)
	{
		// number of items in the list
		int ndata = readInt();
		int i1, i2;
		i1 = i2 = readInt();
		for (int i = 1; i < ndata; i++)
		{
			i2 = readInt();
			if (i2 > 0 && i1 >= 0)
			{
				if (i1 > 0)
				{
					Dof tempVar(ndim * (i1 - 1) + dir, sValue);
					it->add(&tempVar);
				}
				i1 = i2;
			}
			else if (i2 < 0)
			{
				for (int j = i1; j <= (-i2); j++)
				{
					Dof tempVar2(ndim * (j - 1) + dir, sValue);
					it->add(&tempVar2);
				}
				i1 = 0;
				i2 = 0;
			}
		}
		if (i2 > 0)
		{
			Dof tempVar3(ndim * (i2 - 1) + dir, sValue);
			it->add(&tempVar3);
		}
		return it;
	}

	FeScanner::~FeScanner()
	{
		es->close();
	}
}
