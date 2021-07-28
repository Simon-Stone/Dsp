#pragma once

#include <complex>
#include <vector>

#include "utilities.h"
#include "window.h"

namespace dsp
{
	/// @brief Fast Fourier Transformations
	namespace fft
	{
		/// @brief Normalization mode for the various transforms: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		enum class NormalizationMode { backward, ortho, forward };

		/// @brief Backend choices for performing the actual transformations
		enum class backend { automatic, simple, fftw };
		
		/// @brief Frequency range to return (for symmetric transforms)
		enum class FrequencyRange 
		{ 
			centered,  // Return two-sided, centered transform where the negative frequencies are returned first, followed by the positive frequencies
			twosided,  // Return two-sided transform where the positive frequencies are returned first, followed by the negative frequencies
			onesided   // Return one-sided transform where only the positive frequencies are returned (does not conserve the total power!)
		};

		/// @brief Compute the 1-D discrete Fourier Transform.
		///
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) with the efficient Fast Fourier Transform(FFT) algorithm for complex input signals.
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param x Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> cfft(const std::vector<std::complex<T>>& x, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic);

		/// @brief Compute the 1-D inverse discrete Fourier Transform.
		///
		///	This function computes the inverse of the 1-D n-point discrete Fourier transform computed by fft. In other words, ifft(fft(x)) == x to within numerical accuracy.
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param X Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> icfft(const std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic);

		/// @brief Compute the 1-D discrete Fourier Transform for real input.
		/// 
		/// This function computes the 1-D n-point discrete Fourier Transform (DFT) of a real-valued vector by means of an efficient algorithm called the Fast Fourier Transform (FFT).
		/// @tparam T Data type of the real values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param x Real input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.		
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The forward-transformed truncated or zero-padded input.
		template<class T>
		std::vector<std::complex<T>> rfft(const std::vector<T>& x, unsigned n = 0, FrequencyRange range = FrequencyRange::twosided, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic);

		/// @brief Alias for rfft.
		template<class T>
		std::vector<std::complex<T>> fft(const std::vector<T>& x, unsigned n = 0, FrequencyRange range = FrequencyRange::twosided, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic)
		{
			return rfft(x, n, range, mode, backend);
		}
		
		/// @brief Computes the inverse of rfft.
		///
		/// This function computes the inverse of the 1-D n-point discrete Fourier Transform of real input computed by rfft. In other words, irfft(rfft(x), x.size()) == x to within numerical accuracy. 
		/// @tparam T Data type of the complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param X Complex input
		/// @param n Length of the transformed output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is 0 (default), the length of the input is used.
		/// @param mode The normalization mode: "backward" means normalization by n on the inverse transformation only, "forward" means on the forward transformation only, and "ortho" means divide by sqrt(n) in both directions.
		/// @param backend Can be automatic, simple, or fftw. 'simple' is a low-level straight-forward implementation of the complex FFT and 'fftw' uses the FFTW library. 'simple' is best for a small number fo samples due to the overhead of the FFTW planning stage. For longer inputs, FFTW becomes significantly faster. 'automatic' therefore chooses the 'simple' implementation for input lengths of less than 100 000 samples and 'fftw' for longer inputs.
		/// @return The backward-transformed truncated or zero-padded input.
		template<class T>
		std::vector<T> irfft(const std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic);

		/// @brief Alias for irfft()
		template<class T>
		std::vector<T> ifft(const std::vector<std::complex<T>>& X, unsigned n = 0, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic)
		{
			return irfft(X, n, mode, backend);
		}

		/// @brief Computes the logarithmic squared-magnitude spectrum
		/// @tparam T Data type of the samples
		/// @param signal the signal to be transformed
		/// @param N_fft length of the FFT
		/// @param relativeCutoff Determines the number of frequency containers to be returned. The highest frequency index is relativeCutoff * N_fft
		/// @return 
		template<class T>
		std::vector<T> logSquaredMagnitudeSpectrum(const std::vector<T>& signal, int N_fft, double relativeCutoff);
		
		template<class T>
		std::vector<std::complex<T>> dft(const std::vector<T>& x, unsigned n, NormalizationMode mode = NormalizationMode::backward, backend backend = backend::automatic);
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

		/// @brief Computes the Short-Time Fourier Transform (STFT).
		///
		/// This function computes the Short-Time Fourier Transform of a real-valued input signal.
		/// @tparam T Data type of the real and complex values. Should be float, double or long double, other types will cause undefined behavior.
		/// @param x Real input
		/// @param frameLength Length of each frame.
		/// @param overlap Number of overlapping samples of consecutive frames.
		/// @param window Type of window to use before applying the Fourier transform to each frame. Note: Only windows with no parameters (except length) are currently supported!
		/// @param fftLength Length of the FFT of each frame. Must be >= frameLength
		/// @return The complex Short-Time Fourier Transform of the input signal.
		template<class T>
		std::vector<std::vector<std::complex<T>>> stft(const std::vector<T>& x, unsigned frameLength, int overlap, window::type window, unsigned fftLength);

		// Discrete Sin and Cosine Transforms (DST and DCT)
		// TODO:
		// dct();
		//idct();
		//dctn();
		//idctn();
		//dst();
		//idst();
		//dstn();
		//idstn();

		// Helper functions
		// TODO:
		//ifftshift();
		//fftfreq();
		//rfftfreq();
		//next_fast_len();

		/// @brief Fast convolution using the FFT
		template<class T>
		std::vector<T> fftconvolution(const std::vector<T>& volume, const std::vector<T>& kernel, convolution_mode mode = convolution_mode::valid);

		/// @brief Rearranges a Fourier transform X by shifting the zero-frequency component to the center of the vector.
		/// @tparam T Data type of the elements
		/// @param X Vector holding the Fourier transform
		/// @return Fourier transform with the zero-frequency component shifted to the center
		template<class T>
		std::vector<T> fftshift(const std::vector<T>& X);

		/// @brief Rearranges a zero-frequency-shifted Fourier transform X back to the original transform output. In other words, ifftshift undoes the result of fftshift.
		/// @tparam T Data type of the elements
		/// @param X Vector holding the Fourier transform
		/// @return Fourier transform with the zero-frequency component at the lowest index
		template<class T>
		std::vector<T> ifftshift(const std::vector<T>& X);

		/// @brief Returns a spectrogram of the passed signal
		/// @tparam T Data type of the signal's samples
		/// @param signal Signal to analyze
		/// @param frameLength Length of each frame
		/// @param overlap_pct Relative overlap between two frames (e.g., 0.75 for 75% overlap)
		/// @param samplingRate Sampling rate in Hz. Can be -1 to use normalized frequencies.
		/// @param relativeCutoff How much of the spectrum to calculate. Default 0.5 to discard mirrored part of the spectrum (assuming real input).
		/// @param windowType Type of the window to use to window each frame.
		/// @param logSquared Return log-squared magnitude if true (default), otherwise return linear magnitude.
		/// @return A vector containing the log-squared-magnitude spectrum of each windowed frame of the signal.
		template<class T>
		std::vector<std::vector<T>> spectrogram(const std::vector<T>& signal, unsigned frameLength, double overlap_pct = 0.5, 
			int samplingRate = -1, double relativeCutoff = 0.5,  window::type windowType = window::type::hamming, bool logSquared = true);
	}
}
