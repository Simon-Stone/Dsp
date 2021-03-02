#pragma once

#include <vector>

namespace dsp::window
{
	/*
	* The free functions implemented in this module are intended to have a
	* similar interface as the corresponding scipy.signal.windows functions
	*
	*/

	/// <summary>
	/// Types of windows
	/// </summary>
	enum class type {
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
	std::vector<T> boxcar(unsigned N, bool sym = true);

	/// <summary>
	/// Return a triangular window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty vector is returned.</param>
	/// <param name="sym">When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).</returns>
	template<class T>
	std::vector<T> triang(unsigned N, bool sym = true);

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
	std::vector<T> general_cosine(unsigned N, const std::vector<T>& a, bool sym = true);

	/// <summary>
	/// Return a generalized Hamming window.
	/// </summary>	
	template<class T>
	std::vector<T> general_hamming(unsigned N, double alpha, bool sym = true);

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
	std::vector<T> blackman(unsigned N, bool sym = true);

	/// <summary>
	///Return a Hamming window.
	///
	///The Hamming window is a taper formed by using a raised cosine with non - zero endpoints, optimized to minimize the nearest side lobe.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> hamming(unsigned N, bool sym = true);

	/// <summary>
	/// Return a Hann window.
	/// The Hann window is a taper formed by using a raised cosine or sine - squared with ends that touch zero.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> hann(unsigned N, bool sym = true);

	/// <summary>
	/// Return a Bartlett window.
	///
	/// The Bartlett window is very similar to a triangular window, except that the end points are at zero. It is often used in signal processing for tapering a signal, without generating too much ripple in the frequency domain.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The triangular window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> bartlett(unsigned N, bool sym = true);

	/// <summary>
	/// Return a flat top window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> flattop(unsigned N, bool sym = true);

	/// <summary>
	/// Return a Parzen window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> parzen(unsigned N, bool sym = true);

	/// <summary>
	/// Return a Bohman window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> bohman(unsigned N, bool sym = true);

	/// <summary>
	/// Return a minimum 4-term Blackman-Harris window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> blackmanharris(unsigned N, bool sym = true);

	/// <summary>
	/// Return a minimum 4-term Blackman-Harris window according to Nuttall.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> nuttall(unsigned N, bool sym = true);

	/// <summary>
	/// Return a modified Bartlett-Hann window.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> barthann(unsigned N, bool sym = true);

	/// <summary>
	/// Return a Kaiser window.
	///
	/// The Kaiser window is a taper formed by using a Bessel function.
	/// </summary>
	/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="beta">Shape parameter, determines trade-off between main-lobe width and side lobe level. As beta gets large, the window narrows.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> kaiser(unsigned N, double beta, bool sym = true);

	/// <summary>
	/// Return a Gaussian window.
	/// </summary>
/// <typeparam name="T">Type of returned values.</typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="std">The standard deviation, sigma.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is True).</returns>
	template<class T>
	std::vector<T> gaussian(unsigned N, double std, bool sym = true);

	/// <summary>
	/// Return a window with a generalized Gaussian shape.
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="p">Shape parameter. p = 1 is identical to gaussian, p = 0.5 is the same shape as the Laplace distribution.</param>
	/// <param name="sig">The standard deviation, sigma.</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value normalized to 1 (though the value 1 does not appear if M is even and sym is True).</returns>
	template<class T>
	std::vector<T> general_gaussian(unsigned N, double p, double sig, bool sym = true);


	//TODO: dpss window (https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.dpss.html#scipy.signal.windows.dpss)

	/// <summary>
	/// Return a Dolph-Chebyshev window.
	/// </summary>
	/// <typeparam name="T"></typeparam>
	/// <param name="N">Number of points in the output window. If zero or less, an empty array is returned.</param>
	/// <param name="at">Attenuation (in dB).</param>
	/// <param name="sym">When True (default), generates a symmetric window, for use in filter design. When False, generates a periodic window, for use in spectral analysis.</param>
	/// <returns>The window, with the maximum value always normalized to 1</returns>
	template<class T>
	std::vector<T> chebwin(unsigned N, double at, bool sym = true);

}

