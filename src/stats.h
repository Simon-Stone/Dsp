#pragma once
#include <numeric>
#include <vector>

#include "Signal.h"

/// @brief Functions for statistical calculations on std::vector and dsp::Signal
namespace dsp
{
	/// @brief Returns the mean value of the passed range
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @tparam InputIt Iterator type.
	/// @param begin Start of the range.
	/// @param end End of the range.
	/// @return Mean value of the range.
	template<class T, class InputIt>
	T mean(InputIt begin, InputIt end)
	{
		auto sum = std::accumulate(begin, end, T(0));
		return sum / std::distance(begin, end);
	}
	
	/// @brief Returns the mean value of a vector
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param x Vector to calculate the mean of.
	/// @return Mean value of the vector.
	template<class T>
	T mean(const std::vector<T>& x)
	{
		return mean<T>(x.begin(), x.end());
	}

	/// @brief Returns the mean value of a signal
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types may cause undefined behavior.
	/// @param x Signal to calculate the mean of.
	/// @return Mean value of the signal.
	template<class T>
	T mean(const Signal<T>& x)
	{
		return mean<T>(x.begin(), x.end());
	}

	/// @brief Returns the median of a range
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types may cause undefined behavior.
	/// @tparam InputIt Iterator type
	/// @param begin Start of the range
	/// @param end End of the range
	/// @return The median value of the range
	template<class T, class InputIt>
	T median(InputIt begin, InputIt end)
	{
		std::vector<T> tmp(std::distance(begin, end) / 2 + 1);
		std::partial_sort_copy(begin, end, tmp.begin(), tmp.end());
		if (std::distance(begin, end) % 2)
		{
			return tmp.back();
		}
		return static_cast<T>((tmp[tmp.size() - 1] + tmp[tmp.size() - 2]) / 2.0);
	}

	/// @brief Returns the median of a vector or Signal.
	/// @tparam T Type of container. Should be std::vector or dsp::Signal. Other containers may cause undefined behavior.
	/// @param x Container to calculate the median of
	/// @return The median of the container
	template<class T>
	auto median(const T& x)
	{
		return median<typename T::value_type>(x.begin(), x.end());
	}

	/// @brief Calculate the variance or standard deviation based on a sample or based on the population
	enum class weight{sample, population};

	/// @brief Calculates the variance of a range
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @tparam InputIt Iterator type.
	/// @param begin Start of the range.
	/// @param end End of the range.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Variance of the samples in the range.
	template<class T, class InputIt>
	T var(InputIt begin, InputIt end, weight w = weight::sample)
	{
		auto mu = mean<T>(begin, end);
		auto sum = std::transform_reduce(begin, end, T(0), std::plus<T>(), [mu](auto x) {return std::pow(x - mu, 2); });
		if(w == weight::sample)
		{
			return sum / (std::distance(begin, end) - 1);
		}
		return sum / static_cast<T>(std::distance(begin, end));		
	}

	/// @brief Calculates the variance of a vector.
	/// @tparam T Type of the values. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param x Vector to calculate the variance of.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Variance of the samples in the vector.
	template<class T>
	T var(const std::vector<T>& x, weight w = weight::sample)
	{
		return var<T>(x.begin(), x.end(), w);
	}
		
	/// @brief Calculates the variance of a Signal.
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param x Signal to calculate the variance of.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Variance of the signal.
	template<class T>
	T var(const dsp::Signal<T>& x, weight w = weight::sample)
	{
		return var<T>(x.begin(), x.end(), w);
	}

	/// @brief Calculates the standard deviation of a range
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @tparam InputIt Iterator type.
	/// @param begin Start of the range.
	/// @param end End of the range.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Standard deviation of the samples in the range.
	template<class T, class IterIt>
	T std(IterIt begin, IterIt end, weight w = weight::sample)
	{
		return std::sqrt(dsp::var<T>(begin, end, w));
	}

	/// @brief Calculates the standard deviation of a vector.
	/// @tparam T Type of the values. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param x Vector to calculate the standard deviation of.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Standard deviation of the samples in the vector.	
	template<class T>
	T std(const std::vector<T>& x, weight w = weight::sample)
	{
		return std<T>(x.begin(), x.end(), w);
	}

	/// @brief Calculates the standard deviation of a Signal.
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param x Signal to calculate the standard deviation of.
	/// @param w Denominator to use in the calculation: If w == weight::sample, the denominator is N-1; if it w == population, the denominator is N. Default: sample.
	/// @return Standard deviation of the signal.
	template<class T>
	T std(const dsp::Signal<T>& x, weight w = weight::sample)
	{
		return std<T>(x.begin(), x.end(), w);
	}

	/// @brief Returns the standardized z-scores of the data in a range
	/// @tparam InputIt Iterator type
	/// @tparam T Data type of the samples. Should be float, double, or long double. Other types may cause undefined behavior.
	/// @param begin Start of the range
	/// @param end End of the range
	/// @return A vector holding the z-scores of the data in the range (i.e., the standardized data).
	template<class T, class InputIt>
	std::vector<T> zscore(InputIt begin, InputIt end)
	{
		std::vector<T> tmp{ begin, end };
		auto mu = dsp::mean<T>(begin, end);
		auto std = dsp::std<T>(begin, end);
		std::transform(tmp.begin(), tmp.end(), tmp.begin(), [mu, std](auto x) {return (x - mu) / std; });
		return tmp;
	}

	/// @brief Returns the standardized z-scores of the data in a vector
	/// @tparam T Data type of the samples. Should be float, double, or long double. Other types may cause undefined behavior.
	/// @param x Vector that holds the un-standardized data.
	/// @return Vector holding the z-scores of the data in x (i.e., the standardized data). 
	template<class T>
	std::vector<T> zscore(const std::vector<T>& x)
	{
		return zscore<T>(x.begin(), x.end());		
	}

	/// @brief Returns the standardized z-scores of the data in a Signal
	/// @tparam T Data type of the samples. Should be float, double, or long double. Other types may cause undefined behavior.
	/// @param x Signal that holds the un-standardized data.
	/// @return Signal holding the z-scores of the data in x (i.e., the standardized data). 
	template<class T>
	dsp::Signal<T> zscore(const dsp::Signal<T>& x)
	{
		return dsp::Signal<T>(x.getSamplingRate_Hz(), zscore<T>(x.begin(), x.end()));
	}
}