#pragma once
#include <vector>

namespace dsp
{
	/// @brief Special functions
	namespace special
	{
#ifndef ZERO_DEPENDENCIES
		/// @brief Modified Bessel function of order 0.
		///
		/// This function is a wrapper for the boost routine i0.
		/// @tparam T Type of values
		/// @param x Argument
		/// @return Value of the modified Bessel function of order 0 at x
		template<class T>
		T i0(T x);
		
		/// @brief Modified Bessel function of order 0.
		///
		/// This function is a wrapper for the boost routine i0.
		/// @tparam T Type of values
		/// @param x Vector of arguments
		/// @return Value of the modified Bessel function of order 0 at x
		template<class T>
		std::vector<T> i0(std::vector<T> x);
#endif
	}
}
