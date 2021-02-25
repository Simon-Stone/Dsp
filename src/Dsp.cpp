#include "Dsp.h"
#include <cmath>

/// <summary>
/// Return the next-largest power of two k so that n <= 2^k
/// </summary>
/// <param name="n">Number of which find the next-largest power of two.</param>
/// <returns>The next-largest power of two so that n <= 2^k </returns>
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
		sineSignal.push_back(static_cast<T>(std::sin(2.0 * Dsp::pi * frequency_Hz * k * 1.0 / samplingRate_Hz)));
	}
	
	return sineSignal;
}

template Dsp::Signal<float> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<short> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<int> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);


