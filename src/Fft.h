#pragma once
#include "Dsp.h"
#include "Signal.h"

#include <cmath>

namespace Dsp
{
	/*
	 * The free functions implemented in this module are intended to have a
	 * similar interface as the corresponding scipy.fft functions
	 *
	 */

	 // Fast Fourier Transforms (FFT)
	enum class NormalizationMode { backward, ortho, forward };

	/// <summary>
	/// The Fast Fourier Transform for complex-valued signals
	/// </summary>
	/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
	/// <param name="x">The complex-valued original signal.</param>
	/// <param name="n">The length of the FFT. Will be expanded to the next power of two, if it is not already a power of two.</param>
	/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
	/// <returns>A complex-valued signal containing the complex Fourier spectrum.</returns>
	template<class T>
	auto fft(Signal<std::complex<T>>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward)
	{
		// Returned signal is complex number of the same type
		Signal<std::complex<T>> X(x);
		// Default length is the length of the signal
		n = n == 0 ? static_cast<unsigned>(x.size()) : n;
		// Length should be power of two for speed
		const auto exponent = nextpow2(n);
		n = 1 << exponent;
		// Resize signal (will be truncated or zero-padded if necessary)
		X.resize(n, { 0, 0 });

		auto nm1 = n - 1;
		auto nd2 = n / 2;

		auto j = nd2;
		for (unsigned i = 1; i <= n - 2; ++i)
		{
			if (i < j)
			{
				std::swap(X[j], X[i]);
			}
			auto k = nd2;

			while (k <= j)
			{
				j -= k;
				k /= 2;
			}
			j += k;
		}

		for (unsigned l = 1; l <= exponent; ++l)
		{
			const unsigned le{ static_cast<unsigned>(1) << l };
			const unsigned le2{ le / 2 };
			std::complex<T> u{ 1.0, 0.0 };
			std::complex<T> s{ static_cast<T>(std::cos(pi / static_cast<T>(le2))), static_cast<T>(-std::sin(pi / static_cast<T>(le2))) };

			for (j = 1; j <= le2; ++j)
			{
				const auto jm1 = j - 1;
				std::complex<T> t;
				for (unsigned i = jm1; i <= nm1; i += le)
				{
					auto ip = i + le2;
					t = { X[ip].real() * u.real() - X[ip].imag() * u.imag(), X[ip].real() * u.imag() + X[ip].imag() * u.real() };
					X[ip] = X[i] - t;
					X[i] += t;
				}

				t.real(u.real());
				u = { t.real() * s.real() - u.imag() * s.imag(), t.real() * s.imag() + u.imag() * s.real() };
			}
		}
		switch (mode)
		{
		case NormalizationMode::backward:
			// No normalization on the forward transform
			break;
		case NormalizationMode::ortho:
			X /= static_cast<T>(sqrt(n));
			break;
		case NormalizationMode::forward:
			X /= static_cast<T>(n);
			break;
		}
		return X;
	}

	/// <summary>
	/// The Inverse Fast Fourier Transform for complex-valued signals (You probably dont' want this, but irfft()!)
	/// </summary>
	/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
	/// <param name="x">The complex-valued transformed signal.</param>
	/// <param name="n">The length of the IFFT. Will be expanded to the next power of two, if it is not already a power of two.</param>
	/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
	/// <returns>The original complex-valued signal.</returns>
	template<class T>
	auto ifft(Signal<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward)
	{
		// Returned signal is complex number of the same type
		Signal<std::complex<T>> x(X);
		// Default length is the length of the signal
		n = n == 0 ? static_cast<unsigned>(X.size()) : n;
		// Length should be power of two for speed
		const auto exponent = nextpow2(n);
		n = 1 << exponent;
		// Resize signal (will be truncated or zero-padded if necessary)
		x.resize(n, { 0, 0 });

		x = conj(x);

		if (mode == NormalizationMode::backward) { mode = NormalizationMode::forward; }
		x = fft(x, n, mode);

		x = conj(x);

		return x;
	}

	/// <summary>
	/// The Fast Fourier Transform for real-valued signals.
	/// </summary>
	/// <typeparam name="T">Data type of the samples. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
	/// <param name="x">The real-valued original signal.</param>
	/// <param name="n">The length of the FFT. Will be expanded to the next power of two, if it is not already a power of two.</param>
	/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
	/// <returns>A complex-valued signal containing the complex Fourier spectrum.</returns>
	template<class T>
	auto rfft(Signal<T>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward)
	{
		// Returned signal is complex number of the same type
		Signal<std::complex<T>> X;
		X.setSamplingRate_Hz(x.getSamplingRate_Hz());
		X.resize(x.size());
		std::transform(x.begin(), x.end(), X.begin(), [](T re) {return std::complex<T>(re, 0.0); });

		// Default length is the length of the signal
		n = n == 0 ? static_cast<unsigned>(x.size()) : n;
		// Length should be power of two for speed
		const auto exponent = nextpow2(n);
		n = 1 << exponent;
		// Resize signal (will be truncated or zero-padded if necessary)
		X.resize(n, { 0, 0 });

		// Divide in even and odd points
		for (unsigned i = 0; i < n / 2; ++i)
		{
			X[i].real(X[2 * i].real());
			X[i].imag(X[2 * i + 1].real());
		}

		// Calculate complex FFT for n/2 number of points
		X = fft(X, n / 2);
		X.resize(n, 0.0);
		
		auto nm1 = n - 1;
		auto nd2 = n / 2;
		auto n4 = (n / 4) - 1;

		for (unsigned i = 1; i <= n4; ++i)
		{
			unsigned im = nd2 - i;
			unsigned ip2 = i + nd2;
			unsigned ipm = im + nd2;
			X[ip2].real((X[i].imag() + X[im].imag()) / 2);
			X[ipm].real(X[ip2].real());
			X[ip2].imag(-(X[i].real() - X[im].real()) / 2);
			X[ipm].imag(-X[ip2].imag());

			X[i].real((X[i].real() + X[im].real()) / 2);
			X[im].real(X[i].real());
			X[i].imag((X[i].imag() - X[im].imag()) / 2);
			X[im].imag(-X[i].imag());
		}

		X[(n * 3 / 4)] = { X[n / 4].imag(), 0.0 };
		X[nd2] = { X[0].imag(), 0.0 };
		X[n / 4].imag(0.0);
		X[0].imag(0.0);

		auto l = exponent;
		auto le = 1 << l;
		auto le2 = le / 2;
		std::complex<T> u{ 1.0, 0.0 };
		std::complex<T> s{ static_cast<T>(std::cos(Dsp::pi / static_cast<T>(le2))), static_cast<T>(-std::sin(Dsp::pi / static_cast<T>(le2))) };
		std::complex<T> t;
		for (int j = 1; j <= le2; ++j)
		{
			auto jm1 = j - 1;
			for (unsigned i = jm1; i <= nm1; i += le)
			{
				auto ip = i + le2;
				t = { X[ip].real() * u.real() - X[ip].imag() * u.imag(), X[ip].real() * u.imag() + X[ip].imag() * u.real() };
				X[ip] = X[i] - t;
				X[i] += t;
			}

			t.real(u.real());
			u.real(t.real() * s.real() - u.imag() * s.imag());
			u.imag(t.real() * s.imag() + u.imag() * s.real());
		}

		switch (mode)
		{
		case NormalizationMode::backward:
			// Only normalize in backwards direction
			break;
		case NormalizationMode::ortho:
			X /= static_cast<T>(sqrt(n));
			break;
		case NormalizationMode::forward:
			X /= static_cast<T>(n);
			break;
		}

		return X;
	}

	/// <summary>
	/// The Inverse Fast Fourier Transform for real-valued original signals.
	/// </summary>
	/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
	/// <param name="x">The complex-valued transformed signal.</param>
	/// <param name="n">The length of the real IFFT. Will be expanded to the next power of two, if it is not already a power of two.</param>
	/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
	/// <returns>The original real-valued signal.</returns>
	template<class T>
	auto irfft(Signal<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward)
	{
		// Returned signal is a number of the same type
		Signal<T> x(X.getSamplingRate_Hz());
		// Default length is the length of the signal
		const auto n_out = n == 0 ? static_cast<unsigned>(X.size()) : n;
		// Length should be power of two for speed
		const auto exponent = nextpow2(n_out);
		n = 1 << exponent;
		// Resize signal (will be truncated or zero-padded if necessary)
		x.resize(n, 0.0);

		// Force spectrum to be symmetric
		for (auto k = n/2 + 1; k < n; ++k)
		{
			X[k] = std::conj(X[n - k]);
		}

		// Add up real and imaginary part
		for (unsigned k = 0; k < n; ++k)
		{
			x[k] = X[k].real() + X[k].imag();
		}
		auto s = rfft(x, n, NormalizationMode::backward);  // Don't normalize just yet

		// Post-processing
		for (unsigned k = 0; k < n; ++k)
		{
			x[k] = s[k].real() + s[k].imag();
		}

		// Resize to requested output length
		x.resize(n_out, 0.0);

		// Normalize
		switch (mode)
		{
			case NormalizationMode::backward:
				x /= n;
				break;
			case NormalizationMode::ortho:
				x /= static_cast<T>(sqrt(n));
				break;
			case NormalizationMode::forward:
				x /= static_cast<T>(n);
				break;
		}
		return x;
	}

	// TODO:
	//fft2();
	//ifft2();
	//fftn();
	//ifftn();
	//rfft2();
	//irfft2();
	//rfftn();
	//irfftn();
	//hfft();
	//ihfft();
	//hfft2();
	//ihfft2();
	//hfftn();
	//ihfftn();

	// Discrete Sin and Cosine Transforms (DST and DCT)
	// TODO:
	//dct();
	//idct();
	//dctn();
	//idctn();
	//dst();
	//idst();
	//dstn();
	//idstn();

	// Helper functions
	//fftshift();
	//ifftshift();
	//fftfreq();
	//rfftfreq();
	//next_fast_len();
}
