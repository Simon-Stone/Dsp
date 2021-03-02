#pragma once
#include "dsp.h"

namespace dsp
{
	namespace fft
	{
		/*
		 * The free functions implemented in this module are intended to have a
		 * similar interface as the corresponding scipy.fft functions
		 *
		 */

		/// <summary>
		/// Normalization mode for the various transforms: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// </summary>
		enum class NormalizationMode { backward, ortho, forward };

		/// <summary>
		/// Backend choices for performing the actual transformations
		/// </summary>
		enum class backend { automatic, simple, fftw };
		
		/// <summary>
		/// Compute the 1-D discrete Fourier Transform.
		///
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform(FFT) algorithm for complex input signals.
		/// </summary>
		/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
		/// <param name="x">Complex input.</param>
		/// <param name="n">Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.</param>
		/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
		/// <param name="overwrite_x">If true, the contents of x can be destroyed; the default is false.</param>
		/// <returns>The transformed truncated or zero-padded input.</returns>
		template<class T>
		std::vector<std::complex<T>> fft(std::vector<std::complex<T>>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic);

		/// <summary>
		/// Compute the 1-D inverse discrete Fourier Transform.
		///	This function computes the inverse of the 1-D n-point discrete Fourier transform computed by fft. In other words, ifft(fft(x)) == x to within numerical accuracy.
		/// </summary>
		/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
		/// <param name="X">Complex input.</param>
		/// <param name="n">Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.</param>
		/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
		/// <param name="overwrite_X">If true, the contents of X can be destroyed; the default is false.</param>
		/// <returns>The transformed truncated or zero-padded input.</returns>
		template<class T>
		std::vector<std::complex<T>> ifft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic);

		/// <summary>
		/// Compute the 1-D discrete Fourier Transform for real input.
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) of a real-valued vector by means of an efficient algorithm called the Fast Fourier Transform (FFT).
		/// </summary>
		/// <typeparam name="T">Data type of the real values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
		/// <param name="x">Real input</param>
		/// <param name="n">Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.</param>
		/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
		/// <param name="overwrite_x">If true, the contents of x can be destroyed; the default is false.</param>
		/// <returns>The forward-transformed truncated or zero-padded input.</returns>
		template<class T>
		std::vector<std::complex<T>> rfft(std::vector<T>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic);

		/// <summary>
		/// Alias for rfft.
		/// </summary>		
		template<class T>
		std::vector<std::complex<T>> fft(std::vector<T>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic)
		{
			return rfft(x, n, mode, overwrite_x, backend);
		}
		
		/// <summary>
		/// Computes the inverse of rfft.
		///
		/// This function computes the inverse of the 1-D n-point discrete Fourier Transform of real input computed by rfft. In other words, irfft(rfft(x), x.size()) == x to within numerical accuracy. 
		/// </summary>
		/// <typeparam name="T">Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.</typeparam>
		/// <param name="X">Complex input.</param>
		/// <param name="n">Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.</param>
		/// <param name="mode">The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.</param>
		/// <param name="overwrite_X">If true, the contents of x can be destroyed; the default is false.</param>
		/// <returns>The backward-transformed truncated or zero-padded input.</returns>
		template<class T>
		std::vector<T> irfft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic);

		/// <summary>
		/// Alias for irfft()
		/// </summary>
		template<class T>
		std::vector<T> ifft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic)
		{
			return ifft(X, n, mode, overwrite_X, backend);
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
		// TODO:
		//fftshift();
		//ifftshift();
		//fftfreq();
		//rfftfreq();
		//next_fast_len();
	}
}
