#pragma once

#include <string>
#include <stdexcept>
#include "rectangularvectors.h"

namespace util
{


	// JFEM data scanner. Delimiters: blank, =.
	class FeScanner
	{

	private:
		Scanner *es;

		// Constructs FE data scanner.
		// fileIn - name of the file containing data.
	public:
		FeScanner(const std::wstring &fileIn);

		// Returns  true if another token is available.
		virtual bool hasNext();

		 // Returns  true if double is next in input.
		virtual bool hasNextDouble();

		// Gives the next token from this scanner.
		virtual std::wstring next();

		// Gives the next double from this scanner.
		virtual double nextDouble();

		// Reads the next integer.
		// Generates an error if next token is not integer.
		virtual int readInt();

		// Reads the next double.
		// Generates an error if next token is not double.
		virtual double readDouble();

		// Advances the scanner past the current line.
		virtual void nextLine();

		// Moves to line which follows a line with the word.
		virtual void moveAfterLineWithWord(const std::wstring &word);

		// Method reads < nNumbers numbers > and places resulting
		// degrees of freedom in a List data structure.
		// Here numbers is a sequence of the type n1 n2 -n3 ...
		// where n2 -n3 means from n2 to n3 inclusive.
		// it - list iterator.
		// dir - direction (1,2,3).
		// nDim - problem dimension (2/3).
		// sValue - specified value.
		// returns - modified list iterator it.
		virtual ListIterator *readNumberList(ListIterator *it, int dir, int ndim, double sValue);

		// Closes this scanner.
		virtual ~FeScanner();

	};

}
