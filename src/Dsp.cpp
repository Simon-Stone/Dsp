#include "dsp.h"

#include <chrono>
#include <cmath>
#include <functional>

#include "fft.h"

unsigned dsp::nextpow2(unsigned n)
{
	return static_cast<unsigned>(ceil(log2(n)));
}

namespace dsp
{
	template<class T>
	std::vector<T> direct_convolution(const std::vector<T>& in1, const std::vector<T>& in2,
		convolution_mode mode)
	{
		throw std::runtime_error("Direct convolution not yet implemented!");
		return {};
	}
}

template <class T>
std::vector<T> dsp::centered(const std::vector<T>& vec, size_t newSize)
{
	auto currentSize = vec.size();
	auto firstIndex = (currentSize - newSize) / 2;
	auto lastIndex = firstIndex + newSize;
	return { vec.begin() + firstIndex, vec.begin() + lastIndex };
}

template <class T>
std::pair<dsp::convolution_method, std::map<dsp::convolution_method, double>> dsp::choose_conv_method(const std::vector<T>& in1,
	const std::vector<T>& in2, convolution_mode mode, bool measure)
{
	std::map<convolution_method, double> times;
	if (measure)
	{
		// TODO
		throw std::runtime_error("Not implemented yet!");
	}

	// TODO: Calculate if FFT-based method is actually faster than direct method
	return { convolution_method::fft, times };
}

template <class T>
std::vector<T> dsp::correlate(const std::vector<T>& in1, const std::vector<T>& in2,
	convolution_mode mode, convolution_method method)
{
	
	return {};
}

template <class T>
std::vector<T> dsp::convolve(const std::vector<T>& in1, const std::vector<T>& in2,
	convolution_mode mode, convolution_method method)
{
	auto volume = in1;
	auto kernel = in2;

	if (method == convolution_method::automatic)
	{
		auto [chosen_method, timing] = choose_conv_method(volume, kernel, mode);
		method = chosen_method;
	}

	switch (method)
	{
	case convolution_method::direct:
		return direct_convolution(volume, kernel, mode);
	case convolution_method::fft:
		return fft::fftconvolution(volume, kernel, mode);
	default:
		throw std::runtime_error("Unknown convolution method!");
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
template dsp::Signal<float> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<double> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long double> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<short> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<int> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);
template dsp::Signal<long> dsp::sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude, double phase);

template std::vector<float> dsp::centered(const std::vector<float>& vec, size_t newSize);
template std::vector<double> dsp::centered(const std::vector<double>& vec, size_t newSize);
template std::vector<long double> dsp::centered(const std::vector<long double>& vec, size_t newSize);

template std::vector<float> dsp::convolve(const std::vector<float>& in1, const std::vector<float>& in2,
	convolution_mode mode, convolution_method method);
template std::vector<double> dsp::convolve(const std::vector<double>& in1, const std::vector<double>& in2,
	convolution_mode mode, convolution_method method);
template std::vector<long double> dsp::convolve(const std::vector<long double>& in1, const std::vector<long double>& in2,
	convolution_mode mode, convolution_method method);