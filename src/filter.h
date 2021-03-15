#pragma once

#include <vector>

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
}