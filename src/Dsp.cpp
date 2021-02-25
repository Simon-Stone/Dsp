#include "Dsp.h"
#include <cmath>
#include <functional>

unsigned Dsp::nextpow2(unsigned n)
{
	return static_cast<unsigned>(ceil(log2(n)));
}

template <class T>
Dsp::Signal<T> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase)
{
	Dsp::Signal<T> sineSignal(samplingRate_Hz);
	const auto numSamples = length_s * samplingRate_Hz;
	for (unsigned k = 0; k < numSamples; ++k)
	{
		sineSignal.push_back(static_cast<T>(amplitude * std::sin(2.0 * Dsp::pi * frequency_Hz * k * 1.0 / samplingRate_Hz + phase)));
	}
	
	return sineSignal;
}


// Explicit template instantiation
template Dsp::Signal<float> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<short> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<int> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
