#pragma once
#include "Signal.h"

namespace Dsp
{
	// Constants
	const double pi = 3.14159265358979311600;
	
	unsigned nextpow2(unsigned n);
	template <class T>
	Signal<T> sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude = 1.0, double phase = 0.0);
}
