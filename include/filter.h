#pragma once

#include <vector>

#include "Signal.h"

/// @brief Functions to apply and design digital filters
namespace dsp::filter
{
	/// @brief Filters the input data x using a rational transfer function defined by the numerator and denominator coefficients b and a.
	/// @tparam T Data type of the samples
	/// @param b Numerator coefficients
	/// @param a Denominator coefficients
	/// @param x Input signal
	/// @return Filtered signal
	template<class T>
	std::vector<T> filter(std::vector<T> b, std::vector<T> a, 
		const std::vector<T>& x);
	
	/// @brief Returns linear prediction filter coefficients
	/// @tparam T Data type of the samples
	/// @param x Input vector
	/// @param N Order of the linear predictor
	/// @return Linear predictor coefficients
	template<class T>
	std::vector<T> lpc(const std::vector<T>& x, unsigned N);

	/// @brief Perform a median filter on a vector.
	/// @tparam T Data type of the samples
	/// @param x Input vector
	/// @param kernel_size Size of the median filter window. Should be odd! Default: 3
	/// @return Vector containing the median filtered samples (same size as input).
	template<class T>
	std::vector<T> medianfilter(const std::vector<T>& x, size_t kernel_size = 3);

	/// @brief Perform a median filter on a Signal.
	/// @tparam T Data type of the samples
	/// @param x Input Signal
	/// @param kernel_size Size of the median filter window. Should be odd! Default: 3
	/// @return Signal containing the median filtered samples (same size as input).
	template<class T>
	dsp::Signal<T> medianfilter(const Signal<T>& x, typename Signal<T>::size_type kernel_size = 3);
}