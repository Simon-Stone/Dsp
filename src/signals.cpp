#pragma once
#include "signals.h"

#include "dsp.h"

template <class T>
dsp::Signal<T> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase)
{
	Signal<T> sineSignal(samplingRate_Hz);
	const auto numSamples = length_s * samplingRate_Hz;
	for (unsigned k = 0; k < numSamples; ++k)
	{
		sineSignal.push_back(static_cast<T>(amplitude * std::sin(2.0 * dsp::pi * frequency_Hz * k * 1.0 / samplingRate_Hz + phase)));
	}

	return sineSignal;
}

template<class T>
dsp::Signal<T> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase)
{
	return sin<T>(frequency_Hz, length_s, samplingRate_Hz, amplitude, dsp::pi / 2 + phase);
}


// Explicit template instantiation
template dsp::Signal<float> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<double> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long double> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<short> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<int> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long> dsp::signals::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);

template dsp::Signal<float> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<double> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long double> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<short> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<int> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long> dsp::signals::cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
