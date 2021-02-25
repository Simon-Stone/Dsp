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

/// <summary>
/// Return a sampled sinusoid signal
/// </summary>
/// <typeparam name="T">Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.</typeparam>
/// <param name="frequency_Hz">The frequency of the sinusoid in Hertz</param>
/// <param name="length_s">The length of the signal in seconds</param>
/// <param name="samplingRate_Hz">The sampling rate in Hertz</param>
/// <param name="amplitude">Amplitude of the sinusoid</param>
/// <param name="phase">Starting phase angle of the sinusoid in rad</param>
/// <returns>A signal containing the specified sinusoid</returns>
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

template Dsp::Signal<float> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long double> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<short> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<int> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template Dsp::Signal<long> Dsp::sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude, double phase);


