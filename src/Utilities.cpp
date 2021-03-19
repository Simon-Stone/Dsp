#include "utilities.h"

#include <chrono>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "fft.h"
#include "Signal.h"

namespace dsp
{
	template<class T>
	std::vector<T> direct_convolution(const std::vector<T>& in1, const std::vector<T>& in2,
		convolution_mode mode)
	{
		throw std::runtime_error("Direct convolution not yet implemented!");
		return {};
	}

	/// @brief Reverse and conjugate a vector
	template<class T>
	std::vector<T> _reverse_and_conj(const std::vector<T>& vec)
	{
		std::vector<T> reversed_copy(vec.size());
		std::reverse_copy(vec.begin(), vec.end(), reversed_copy.begin());
		return dsp::conj(reversed_copy);
	}

	unsigned maxNumFrames(unsigned signalLength, unsigned frameLength, unsigned overlap)
	{
		const auto frameStride = frameLength - overlap;
		return static_cast<unsigned>(floor(static_cast<double>(signalLength) / frameStride + 1));
	}
}

template <class T>
T dsp::calculateEnergy(typename std::vector<T>::iterator start,
	typename std::vector<T>::iterator end)
{
	return static_cast<T>(std::inner_product(start, end, start, T{ 0.0 }));
}

template <class T>
T dsp::calculateMeanPower(typename std::vector<T>::iterator start,
	typename std::vector<T>::iterator end)
{
	auto energy = calculateEnergy<T>(start, end);
	auto numSamples = std::distance(start, end);
	return static_cast<T>(energy / numSamples);
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
std::vector<T> dsp::correlate(const std::vector<T>& in1, const std::vector<T>& in2,
	convolution_mode mode, convolution_method method)
{
	switch (method)
	{
	case convolution_method::automatic:
		// Fall through to fft method
	case convolution_method::fft:
		return convolve(in1, _reverse_and_conj(in2), mode, method);
	case convolution_method::direct:
		throw std::runtime_error("Direct correlation method not yet implemented!");
	default:
		throw std::runtime_error("Unknown correlation method!");
	}
}

template <class T>
std::vector<std::vector<T>> dsp::signalToFrames(const std::vector<T>& signal, unsigned frameLength, unsigned overlap)
{
	std::vector<std::vector<T>> framedSignal;
	const auto numSamples = signal.size();

	const int frameStride = frameLength - overlap;
	for (size_t startSample = 0; startSample < numSamples; startSample += frameStride)
	{
		const size_t finalSample = std::min(startSample + frameLength, signal.size());
		std::vector<T> frame{ signal.begin() + startSample, signal.begin() + finalSample };
		frame.resize(frameLength, 0);  // Make sure the frame is padded to the frame length with zeros
		framedSignal.push_back(frame);
	}

	return framedSignal;
}

std::pair<unsigned, bool> dsp::window::utilities::extend(unsigned N, bool sym)
{
	if (!sym)
	{
		return { N + 1, true };
	}

	return { N, false };
}

// Explicit template instantiation
template std::vector<float> dsp::centered(const std::vector<float>& vec, size_t newSize);
template std::vector<double> dsp::centered(const std::vector<double>& vec, size_t newSize);
template std::vector<long double> dsp::centered(const std::vector<long double>& vec, size_t newSize);

template std::vector<float> dsp::convolve(const std::vector<float>& in1, const std::vector<float>& in2,
	convolution_mode mode, convolution_method method);
template std::vector<double> dsp::convolve(const std::vector<double>& in1, const std::vector<double>& in2,
	convolution_mode mode, convolution_method method);
template std::vector<long double> dsp::convolve(const std::vector<long double>& in1, const std::vector<long double>& in2,
	convolution_mode mode, convolution_method method);

template std::vector<float> dsp::correlate(const std::vector<float>& in1, const std::vector<float>& in2,
	correlation_mode mode, correlation_method method);
template std::vector<double> dsp::correlate(const std::vector<double>& in1, const std::vector<double>& in2,
	correlation_mode mode, correlation_method method);
template std::vector<long double> dsp::correlate(const std::vector<long double>& in1, const std::vector<long double>& in2,
	correlation_mode mode, correlation_method method);

template float dsp::calculateEnergy(std::vector<float>::iterator start,
	std::vector<float>::iterator end);
template double dsp::calculateEnergy(std::vector<double>::iterator start,
	std::vector<double>::iterator end);
template long double dsp::calculateEnergy(std::vector<long double>::iterator start,
	std::vector<long double>::iterator end);

template float dsp::calculateMeanPower(typename std::vector<float>::iterator start,
	typename std::vector<float>::iterator end);
template double dsp::calculateMeanPower(typename std::vector<double>::iterator start,
	typename std::vector<double>::iterator end);
template long double dsp::calculateMeanPower(typename std::vector<long double>::iterator start,
	typename std::vector<long double>::iterator end);

template std::vector<std::vector<float>> dsp::signalToFrames(const std::vector<float>& signal, unsigned frameLength, unsigned overlap);
template std::vector<std::vector<double>> dsp::signalToFrames(const std::vector<double>& signal, unsigned frameLength, unsigned overlap);
template std::vector<std::vector<long double>> dsp::signalToFrames(const std::vector<long double>& signal, unsigned frameLength, unsigned overlap);