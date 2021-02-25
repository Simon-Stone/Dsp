#pragma once

#include <stdexcept>
#include <vector>

#include "Utilities.h"

namespace Dsp::Window
{
	/*
	* The free functions implemented in this module are intended to have a
	* similar interface as the corresponding scipy.signal.windows functions
	*
	*/

	/// <summary>
	/// Types of windows
	/// </summary>
	enum class window {
		boxcar, triang, blackman, hamming, hann, bartlett, flattop,
		parzen, bohman, blackmanharris, nuttal, barthann, kaiser, gaussian,
		general_gaussian, dpss, chebwin, exponential, tukey, taylor
	};

	/// <summary>
	/// Return a boxcar or rectangular window.
	///
	/// Also known as a rectangular window or Dirichlet window, this is equivalent to no window at all.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the window. If zero or less, an empty vector is returned.</param>
	/// <param name="sym">Whether the window is symmetric. Has no effect for boxcar.</param>
	/// <returns>The window, with the maximum value normalized to 1.</returns>
	template<class T>
	std::vector<T> boxcar(unsigned N, bool sym = true)
	{
		if (N <= 0) return std::vector<T>();

		return std::vector<T>(N, 1);
	}

	/// <summary>
	/// Return a triangular window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty vector is returned.</param>
	/// <param name="sym">When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).</returns>
	template<class T>
	std::vector<T> triang(unsigned N, bool sym = true)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);
		auto n = Dsp::arange(1.0, static_cast<T>((M + 1) / 2 + 1));

		std::vector<T> w;
		w.resize(n.size(), 0.0);
		if (M % 2 == 0)
		{
			std::transform(n.begin(), n.end(),
				w.begin(), 
				[M](auto n) {return (2 * n - 1.0) / M; });
			w.insert(w.end(), w.rbegin(), w.rend());
		}
		else
		{
			std::transform(n.begin(), n.end(),
				w.begin(),
				[M](auto n) {return 2 * n / (M + 1.0); });
			w.insert(w.end(), w.rbegin() + 1, w.rend());
		}

		return Utilities::truncate(w, needs_trunc);
	}

	/// <summary>
	/// Generic weighted sum of cosine terms window
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N"> Number of points in the output window. If zero or less, an empty vector is returned.</param>
	/// <param name="a"> Sequence of weighting coefficients. This uses the convention of being centered on
	/// the origin, so these will typically all be positive numbers, not alternating sign.</param>
	/// <param name="sym">When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.</param>
	/// <returns></returns>
	template<class T>
	std::vector<T> general_cosine(unsigned N, const std::vector<T>& a, bool sym = true)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);

		auto fac = Dsp::linspace<T>(-pi, pi, M);
		std::vector<T> w(M, 0);
		for (size_t k = 0; k < a.size(); ++k)
		{
			auto ak = a[k];
			std::transform(w.begin(), w.end(),
				fac.begin(),
				w.begin(),
				[k, ak](auto w, auto fac){return w + ak * std::cos(k * fac); });
		}
		return Utilities::truncate(w, needs_trunc);
	}

	
	/// <summary>
	/// Return a Blackman window.
	///
	/// The Blackman window is a taper formed by using the first three terms of a summation of cosines. It was
	/// designed to have close to the minimal leakage possible. It is close to optimal, only slightly worse
	/// than a Kaiser window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).</returns>
	template<class T>
	std::vector<T> blackman(unsigned N, bool sym = true)
	{
		return general_cosine<T>(N, { 0.42, 0.50, 0.08 }, sym);
	}
	
	template<class T, class ... Types>
	std::vector<T> getWindow(std::tuple<Types ...> window, unsigned N, bool fftbins)
	{
		switch (std::get<0>(window))
		{
		case window::boxcar:
			return boxcar<T>(N, !fftbins);
		default:
			throw std::logic_error("Window type is not implemented yet!");
		}
	}

}

