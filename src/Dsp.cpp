#include "dsp.h"
#include <cmath>
#include <functional>

unsigned dsp::nextpow2(unsigned n)
{
	return static_cast<unsigned>(ceil(log2(n)));
}

template <class T>
dsp::Signal<T> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase)
{
	dsp::Signal<T> sineSignal(samplingRate_Hz);
	const auto numSamples = length_s * samplingRate_Hz;
	for (unsigned k = 0; k < numSamples; ++k)
	{
		sineSignal.push_back(static_cast<T>(amplitude * std::sin(2.0 * dsp::pi * frequency_Hz * k * 1.0 / samplingRate_Hz + phase)));
	}
	
	return sineSignal;
}


// Explicit template instantiation
template dsp::Signal<float> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<double> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long double> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<short> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<int> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long> dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
