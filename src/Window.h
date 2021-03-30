#pragma once

#include <vector>

/// @brief Window functions
namespace dsp::window
{
	/// @brief Types of windows (not all implemented yet)
	enum class type {
		boxcar, triang, blackman, hamming, hann, bartlett, flattop,
		parzen, bohman, blackmanharris, nuttal, barthann, kaiser, gaussian,
		general_gaussian, dpss, chebwin, exponential, tukey, taylor
	};

	template<class T>
	std::vector<T> get_window(type type, unsigned N, bool sym = true, const std::vector<T>& parameters = {});

	/// @brief Return a boxcar or rectangular window
	///
	/// Also known as a rectangular window or Dirichlet window, this is equivalent to no window at all.
	/// @tparam T Type of returned values.
	/// @param N Number of points in the window. If zero or less, an empty vector is returned.
	/// @param sym Whether the window is symmetric. Has no effect for boxcar.
	/// @return The window, with the maximum value normalized to 1.
	template<class T>
	std::vector<T> boxcar(unsigned N, bool sym = true);

	/// @brief Return a triangular window
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty vector is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> triang(unsigned N, bool sym = true);

	/// @brief Generic weighted sum of cosine terms window
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty vector is returned.
	/// @param a Sequence of weighting coefficients. This uses the convention of being centered on
	/// the origin, so these will typically all be positive numbers, not alternating sign.
	/// @param sym >When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> general_cosine(unsigned N, const std::vector<T>& a, bool sym = true);

	/// @brief Return a generalized Hamming window.
	template<class T>
	std::vector<T> general_hamming(unsigned N, double alpha, bool sym = true);

	/// @brief Return a Blackman window.
	///
	/// The Blackman window is a taper formed by using the first three terms of a summation of cosines. It was
	/// designed to have close to the minimal leakage possible. It is close to optimal, only slightly worse
	/// than a Kaiser window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false,
	/// generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> blackman(unsigned N, bool sym = true);

	/// @brief Return a Hamming window.
	///
	///The Hamming window is a taper formed by using a raised cosine with non - zero endpoints, optimized to minimize the nearest side lobe.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> hamming(unsigned N, bool sym = true);

	/// @brief Return a Hann window.
	///
	/// The Hann window is a taper formed by using a raised cosine or sine - squared with ends that touch zero.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> hann(unsigned N, bool sym = true);

	/// @brief Return a Bartlett window.
	///
	/// The Bartlett window is very similar to a triangular window, except that the end points are at zero. It is often used in signal processing for tapering a signal, without generating too much ripple in the frequency domain.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The triangular window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> bartlett(unsigned N, bool sym = true);

	/// @brief Return a flat top window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> flattop(unsigned N, bool sym = true);

	/// @brief Return a Parzen window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> parzen(unsigned N, bool sym = true);

	/// @brief Return a Bohman window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> bohman(unsigned N, bool sym = true);

	/// @brief Return a minimum 4-term Blackman-Harris window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> blackmanharris(unsigned N, bool sym = true);

	/// @brief Return a minimum 4-term Blackman-Harris window according to Nuttall.
	/// @tparam T Type of returned values.
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> nuttall(unsigned N, bool sym = true);

	/// @brief Return a modified Bartlett-Hann window.
	/// @tparam T Type of returned values.
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> barthann(unsigned N, bool sym = true);

#ifndef ZERO_DEPENDENCIES
	/// @brief Return a Kaiser window.
	/// 
	/// The Kaiser window is a taper formed by using a Bessel function.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param beta Shape parameter, determines trade-off between main-lobe width and side lobe level. As beta gets large, the window narrows.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> kaiser(unsigned N, double beta, bool sym = true);
#endif
	/// @brief Return a Gaussian window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param std The standard deviation, sigma.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the first and last samples equal to zero and the maximum value normalized to 1 (though the value 1 does not appear if N is even and sym is true).
	template<class T>
	std::vector<T> gaussian(unsigned N, double std, bool sym = true);

	/// @brief Return a window with a generalized Gaussian shape.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param p Shape parameter. p = 1 is identical to gaussian, p = 0.5 is the same shape as the Laplace distribution.
	/// @param sig The standard deviation, sigma
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if M is even and sym is true).
	template<class T>
	std::vector<T> general_gaussian(unsigned N, double p, double sig, bool sym = true);

	//TODO: dpss window (https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.dpss.html#scipy.signal.windows.dpss)

	/// @brief Return a Dolph-Chebyshev window
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param at Attenuation (in dB).
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value always normalized to 1
	template<class T>
	std::vector<T> chebwin(unsigned N, T at, bool sym = true);

	/// @brief Return an exponential (or Poisson) window.
	/// @tparam T Type of returned values
	/// @param center Parameter defining the center location of the window function. The default value if not given is center = (M-1) / 2. This parameter must take its default value for symmetric windows.
	/// @param N Number of points in the output window. If zero or less, an empty array is returned.
	/// @param tau Parameter defining the decay. For center = 0 use tau = -(M-1) / ln(x) if x is the fraction of the window remaining at the end.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if M is even and sym is true).
	template<class T>
	std::vector<T> exponential(T center, unsigned N, T tau = 1.0, bool sym = true);

	/// @brief Alias for exponential() that provides the default value for center (which depends on N)
	template<class T>
	std::vector<T> exponential(unsigned N, T tau = 1.0, bool sym = true)
	{
		return exponential<T>((N - 1) / 2, N, tau, sym);
	}

	/// @brief Return a Tukey window, also known as a tapered cosine window.
	/// @tparam T Type of returned values
	/// @param N Number of points in the output window. If zero or less, an empty vector is returned.
	/// @param alpha Shape parameter of the Tukey window, representing the fraction of the window inside the cosine tapered region. If zero, the Tukey window is equivalent to a rectangular window. If one, the Tukey window is equivalent to a Hann window.
	/// @param sym When true (default), generates a symmetric window, for use in filter design. When false, generates a periodic window, for use in spectral analysis.
	/// @return The window, with the maximum value normalized to 1 (though the value 1 does not appear if M is even and sym is true).
	template <class T>
	std::vector<T> tukey(unsigned N, T alpha = 0.5, bool sym = true);


	// TODO: taylor window https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.taylor.html#scipy.signal.windows.taylor
}

