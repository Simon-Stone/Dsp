#pragma once
#include "Signal.h"

/// @brief Convenience functions for digital signal processing
namespace dsp
{
	// Constants
	const double pi = 3.14159265358979311600;

	/// @brief Return the next-largest power of two k so that n <= 2^k
	/// @param n Number of which find the next-largest power of two.
	/// @return The next-largest power of two so that n <= 2^k 
	unsigned nextpow2(unsigned n);

	/// @brief Return the elements of a vector that satisfy some condition.
	/// @tparam T 
	/// @tparam Pred 
	/// @param condition A predicate which if true indicates the elements of vec to extract.
	/// @param vec Input vector
	/// @return Vector of values from vec where condition is true.
	template <class T, class Pred>
	std::vector<T> extract(Pred condition, const std::vector<T>& vec)
	{
		std::vector<T> results;

		auto it = std::find_if(vec.begin(), vec.end(), condition);
		while (it != vec.end())
		{
			results.push_back(*it);
			it = std::find_if(std::next(it), vec.end(), condition);
		}
		return results;
	}

	/// @brief Return evenly spaced values within a given interval.
	///
	/// Values are generated within the half-open interval [start, stop) (in other
	/// words, the interval including start but excluding stop). 
	/// @tparam T Type of the elements in vector.
	/// @param start Start of interval. The interval includes this value.
	/// @param stop End of interval. The interval does not include this value.
	/// @param step Spacing between values. For any output out, this is the distance
	/// between two adjacent values, out[i+1] - out[i]. The default step size is 1. 
	/// @return Vector of evenly spaced values
	template<typename T>
	std::vector<T> arange(T start, T stop, T step = 1)
	{
		std::vector<T> values;
		for (T value = start; value < stop; value += step)
		{
			values.push_back(value);
		}
		return values;
	}

	/// @brief Return evenly spaced numbers over a specified interval.
	///
	/// Returns N evenly spaced samples, calculated over the interval[start, stop].
	///
	/// The endpoint of the interval can optionally be excluded.
	/// @tparam T Type of the values in the vector.
	/// @param start The starting value of the sequence.
	/// @param stop The end value of the sequence, unless endpoint is set to False. In that case, the sequence consists of all but the last of num + 1 evenly spaced samples, so that stop is excluded. Note that the step size changes when endpoint is False.
	/// @param num Number of samples to generate. Default is 50.
	/// @param endpoint If True, stop is the last sample. Otherwise, it is not included. Default is True.
	/// @return Vector of evenly spaced numbers over the specified interval
	template <typename T>
	std::vector<T> linspace(T start, T stop, size_t num = 50, bool endpoint = true)
	{
		T h = (stop - start) / static_cast<T>(num - 1);
		std::vector<T> xs(num);

		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = start; x != xs.end(); ++x, val += h)
		{
			*x = val;
		}
		if (!endpoint)
		{
			xs.pop_back();
		}
		return xs;
	}

	/// @brief Return a window of a given length and type.
	/// @tparam T Data type of the samples.
	/// @tparam ...Args Types of the elements in window. First is always dsp::window::type, others are optional (see below).
	/// @param window First element must be the window type. If the window requires more parameters, they must be provided here as further elements.
	/// @param Nx The number of samples in the window.
	/// @param fftbins If true (default), create a “periodic” window, ready to use with ifftshift and be multiplied by the result of an FFT (see also fftfreq). If false, create a “symmetric” window, for use in filter design.
	/// @return Returns a window of length Nx and type window.
	template<class T, class... Args>
	std::vector<T> get_window(std::tuple<Args...> window, unsigned Nx, bool fftbins = true);
	
	/// @brief Return a sampled sinusoid signal
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param frequency_Hz The frequency of the sinusoid in Hertz
	/// @param length_s The length of the signal in seconds
	/// @param samplingRate_Hz The sampling rate in Hertz
	/// @param amplitude Amplitude of the sinusoid
	/// @param phase Starting phase angle of the sinusoid in rad
	/// @return A signal containing the specified sinusoid
	template <class T>
	Signal<T> sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude = 1.0, double phase = 0.0);
}
