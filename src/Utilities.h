#pragma once
#include <utility>
#include <vector>

namespace dsp
{
	/// <summary>
	/// Return evenly spaced values within a given interval.
	///
	/// Values are generated within the half-open interval [start, stop) (in other
	/// words, the interval including start but excluding stop). 
	/// </summary>
	/// <typeparam name="T">Type of the elements in vector.</typeparam>
	/// <param name="start">Start of interval. The interval includes this value.</param>
	/// <param name="stop">nd of interval. The interval does not include this value.</param>
	/// <param name="step">Spacing between values. For any output out, this is the distance
	/// between two adjacent values, out[i+1] - out[i]. The default step size is 1. </param>
	/// <returns></returns>
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

	/// <summary>
	/// Return evenly spaced numbers over a specified interval.
	///
	/// Returns N evenly spaced samples, calculated over the interval[start, stop].
	///
	/// The endpoint of the interval can optionally be excluded.
	/// </summary>
	/// <typeparam name="T">Type of the values in the vector.</typeparam>
	/// <param name="start">The starting value of the sequence.</param>
	/// <param name="stop">The end value of the sequence, unless endpoint is set to False. In that case, the sequence consists of all but the last of num + 1 evenly spaced samples, so that stop is excluded. Note that the step size changes when endpoint is False.</param>
	/// <param name="num">Number of samples to generate. Default is 50.</param>
	/// <param name="endpoint">If True, stop is the last sample. Otherwise, it is not included. Default is True.</param>
	/// <returns></returns>
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
		if(!endpoint)
		{
			xs.pop_back();
		}
		return xs;
	}
	
	namespace window
	{
		namespace utilities
		{
			/// <summary>
			/// Extend window length by 1 sample if needed for DFT-even symmetry
			/// </summary>
			std::pair<unsigned, bool> extend(unsigned N, bool sym);

			/// <summary>
			/// Truncate window by 1 sample if needed for DFT-even symmetry
			/// </summary>
			template<class T>
			std::vector<T> truncate(std::vector<T>& w, bool needed)
			{
				if (needed)
				{
					w.pop_back();
					return w;
				}
				return w;					
			}			
		}
	}	
}
