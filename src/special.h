#pragma once
#include <vector>

namespace dsp
{
	namespace special
	{
		/// <summary>
		/// Modified Bessel function of order 0.
		///
		/// This function is a wrapper for the cephes routine i0.
		/// </summary>
		/// /// <typeparam name="T">Type of values.</typeparam>
		/// <param name="x">Argument</param>
		/// <returns>Value of the modified Bessel function of order 0 at x.</returns>
		template<class T>
		T i0(T x);
		
		/// <summary>
		/// Modified Bessel function of order 0.
		/// </summary>
		/// <typeparam name="T">Type of values.</typeparam>
		/// <param name="x">Vector of arguments</param>
		/// <returns>Vector of values of the modified Bessel function of order 0 at x.</returns>
		template<class T>
		std::vector<T> i0(std::vector<T> x);
	}
}
