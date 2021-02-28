#pragma once
#include "Signal.h"

namespace dsp
{
	// Constants
	const double pi = 3.14159265358979311600;

	/// <summary>
	/// Return the next-largest power of two k so that n <= 2^k
	/// </summary>
	/// <param name="n">Number of which find the next-largest power of two.</param>
	/// <returns>The next-largest power of two so that n <= 2^k </returns>
	unsigned nextpow2(unsigned n);

	/// <summary>
	/// Return the elements of a vector that satisfy some condition.
	/// </summary>
	/// <typeparam name="UnaryPred"></typeparam>
	/// <typeparam name="T"></typeparam>
	/// <param name="condition">A predicate which if true indicates the elements of vec to extract.</param>
	/// <param name="vec">Input vector</param>
	/// <returns>Vector of values from vec where condition is true.</returns>
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
	Signal<T> sin(unsigned frequency_Hz, unsigned length_s, unsigned samplingRate_Hz, double amplitude = 1.0, double phase = 0.0);
}
