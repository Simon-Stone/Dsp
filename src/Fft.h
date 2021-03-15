#pragma once

#include <complex>
#include <vector>

#include "utilities.h"


namespace dsp
{
	/// @brief Fast Fourier Transformations
	namespace fft
	{
		/// @brief Normalization mode for the various transforms: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		enum class NormalizationMode { backward, ortho, forward };

		/// @brief Backend choices for performing the actual transformations
		enum class backend { automatic, simple, fftw };
		
		/// @brief Compute the 1-D discrete Fourier Transform.
		///
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform(FFT) algorithm for complex input signals.
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param x Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param overwrite_x If true, the contents of x can be destroyed; the default is false.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> cfft(std::vector<std::complex<T>>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic);

		/// @brief Compute the 1-D inverse discrete Fourier Transform.
		///
		///	This function computes the inverse of the 1-D n-point discrete Fourier transform computed by fft. In other words, ifft(fft(x)) == x to within numerical accuracy.
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param X Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param overwrite_X If true, the contents of X can be destroyed; the default is false.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> icfft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic);

		/// @brief Compute the 1-D discrete Fourier Transform for real input.
		/// 
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) of a real-valued vector by means of an efficient algorithm called the Fast Fourier Transform (FFT).
		/// @tparam T Data type of the real values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param x Real input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param overwrite_x If true, the contents of x can be destroyed; the default is false.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The forward-transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> rfft(std::vector<T>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic);

		/// @brief Alias for rfft.
		template<class T>
		std::vector<std::complex<T>> fft(std::vector<T>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_x = false, backend backend = backend::automatic)
		{
			return rfft(x, n, mode, overwrite_x, backend);
		}
		
		/// @brief Computes the inverse of rfft.
		///
		/// This function computes the inverse of the 1-D n-point discrete Fourier Transform of real input computed by rfft. In other words, irfft(rfft(x), x.size()) == x to within numerical accuracy. 
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param X Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param overwrite_X If true, the contents of x can be destroyed; the default is false.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The backward-transformed truncated or zero-padded input.
		template<class T>
		std::vector<T> irfft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic);

		/// @brief Alias for irfft()
		template<class T>
		std::vector<T> ifft(std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, bool overwrite_X = false, backend backend = backend::automatic)
		{
			return irfft(X, n, mode, overwrite_X, backend);
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

		/// @brief Fast convolution using the FFT
		template<class T>
		std::vector<T> fftconvolution(const std::vector<T>& volume, const std::vector<T>& kernel, convolution_mode mode = convolution_mode::valid);
	}
}
