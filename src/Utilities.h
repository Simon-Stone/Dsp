#pragma once
#include <utility>
#include <vector>

namespace dsp
{
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
			/// @brief Extend window length by 1 sample if needed for DFT-even symmetry		
			std::pair<unsigned, bool> extend(unsigned N, bool sym);

			/// @brief Truncate window by 1 sample if needed for DFT - even symmetry
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
