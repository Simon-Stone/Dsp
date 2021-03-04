#include "dsp.h"

#include <cmath>
#include <functional>

#include "window.h"

unsigned dsp::nextpow2(unsigned n)
{
	return static_cast<unsigned>(ceil(log2(n)));
}

template <class T, class ... Args>
std::vector<T> dsp::get_window(std::tuple<Args...> window, unsigned Nx, bool fftbins)
{
	if (Nx == 0) { return {}; }

	dsp::window::type type = std::get<0>(window);
	switch (type)
	{
	case window::type::boxcar:
		return window::boxcar<T>(Nx, !fftbins);
	case window::type::triang:
		return window::triang<T>(Nx, !fftbins);
	case window::type::blackman:
		return window::blackman<T>(Nx, !fftbins);
	case window::type::hamming:
		return window::hamming<T>(Nx, !fftbins);
	case window::type::hann:
		return window::hann<T>(Nx, !fftbins);
	case window::type::bartlett:
		return window::bartlett<T>(Nx, !fftbins);
	case window::type::flattop:
		return window::flattop<T>(Nx, !fftbins);
	case window::type::parzen:
		return window::parzen<T>(Nx, !fftbins);
	case window::type::bohman:
		return window::bohman<T>(Nx, !fftbins);
	case window::type::blackmanharris:
		return window::blackmanharris<T>(Nx, !fftbins);
	case window::type::nuttal:
		return window::nuttall<T>(Nx, !fftbins);
	case window::type::barthann:
		return window::barthann<T>(Nx, !fftbins);
	case window::type::kaiser:
		return window::barthann<T>(Nx, !fftbins);
	case window::type::gaussian:
		return window::gaussian<T>(Nx, std::get<1>(window), !fftbins);
	case window::type::general_gaussian:
		return window::general_gaussian<T>(Nx, std::get<1>(window), !fftbins);
	case window::type::chebwin:
		return window::chebwin<T>(Nx, std::get<1>(window), !fftbins);
	case window::type::dpss:
	case window::type::exponential:
	case window::type::tukey:
	case window::type::taylor:
		throw std::runtime_error("Window type not yet implemented!");
	default:
		throw std::runtime_error("Unknown window type requested!");
	}

}

template <class T>
dsp::Signal<T> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase)
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
//template std::vector<float> dsp::get_window(std::tuple<window::type, double> window, unsigned Nx, bool fftbins);
template std::vector<double> dsp::get_window(std::tuple<window::type, double> window, unsigned Nx, bool fftbins);
//template std::vector<long double> dsp::get_window(std::tuple<window::type, double> window, unsigned Nx, bool fftbins);
template dsp::Signal<float> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<double> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long double> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<short> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<int> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
