#include "Utilities.h"

std::pair<unsigned, bool> Dsp::Window::Utilities::extend(unsigned N, bool sym)
{
	if (!sym)
	{
		return { N + 1, true };
	}

	return { N, false };
}