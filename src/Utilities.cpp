#include "utilities.h"

std::pair<unsigned, bool> dsp::window::utilities::extend(unsigned N, bool sym)
{
	if (!sym)
	{
		return { N + 1, true };
	}

	return { N, false };
}