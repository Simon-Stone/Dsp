#include "fft.h"

#include <algorithm>
#include <execution>

#ifndef ZERO_DEPENDENCIES
#include <fftw3/fftw3.h>
#endif

#include "filter.h"
#include "Signal.h"
#include "signals.h"


namespace dsp::fft
{
	/// @cond developer-only

	/// @brief Helper function to get the FFT length
	template<class T>
	auto get_fft_length(const std::vector<T>& x, unsigned n)
	{
		if (n == 0)
		{
			n = static_cast<unsigned>(x.size());
		}
		return 2 << (nextpow2(n) - 1);
	}

	template<class T>
	auto resize_fft_input(std::vector<T> x, unsigned n)
	{
		x.resize(n, T());
		return x;
	}

	/* Straight-forward FFT implementations (high-performance for short signals) */
	template<class T>
	std::vector<std::complex<T>> fft_(const std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode)
	{
		if (x.empty()) return {};

		int N = get_fft_length(x, n);

		std::vector<std::complex<T>> X = resize_fft_input(x, N);
		auto* in = reinterpret_cast<T*>(&X[0]);

		int nm1 = N - 1;
		int nd2 = N / 2;
		int j = nd2;
		for (int i = 1; i <= N - 2; ++i)
		{
			if (i < j)
			{
				T t[2] = { in[2 * j], in[2 * j + 1] };
				in[2 * j] = in[2 * i];
				in[2 * j + 1] = in[2 * i + 1];
				in[2 * i] = t[0];
				in[2 * i + 1] = t[1];
			}
			int k = nd2;

			while (k <= j)
			{
				j -= k;
				k /= 2;
			}
			j += k;
		}
		auto exponent = static_cast<unsigned>(std::log2(N));
		for (unsigned l = 1; l <= exponent; ++l)
		{
			auto le = 1 << l;
			auto le2 = le / 2;
			T u[2] = { 1.0, 0.0 };
			T s[2] = { static_cast<T>(std::cos(pi / le2)), static_cast<T>(-1 * std::sin(pi / le2)) };
			T t[2];
			for (auto j = 1; j <= le2; ++j)
			{
				auto jm1 = j - 1;

				for (auto i = jm1; i <= nm1; i += le)
				{
					auto ip = i + le2;
					t[0] = in[2 * ip] * u[0] - in[2 * ip + 1] * u[1];
					t[1] = in[2 * ip] * u[1] + in[2 * ip + 1] * u[0];
					in[2 * ip] = in[2 * i] - t[0];
					in[2 * ip + 1] = in[2 * i + 1] - t[1];
					in[2 * i] += t[0];
					in[2 * i + 1] += t[1];
				}

				t[0] = u[0];
				u[0] = t[0] * s[0] - u[1] * s[1];
				u[1] = t[0] * s[1] + u[1] * s[0];
			}
		}

		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<T>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<T>(N); });
			break;
		}

		return X;
	}

	template<class T>
	std::vector<std::complex<T>> ifft_(const std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode)
	{
		if (x.empty()) return {};

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<T>> X = resize_fft_input(x, N);

		X = dsp::conj(X);

		// Switch the mode because we are using the forward transform in the backward transform
		if (mode == NormalizationMode::backward)
		{
			mode = NormalizationMode::forward;
		}
		else if (mode == NormalizationMode::forward)
		{
			mode = NormalizationMode::backward;
		}
		X = fft_(X, N, mode);

		return dsp::conj(X);
	}

	template<class T>
	std::vector<std::complex<T>> rfft_(const std::vector<T>& x, unsigned n, NormalizationMode mode)
	{
		if (x.empty()) return {};

		unsigned N = get_fft_length(x, n);

		std::vector<T> x_in = resize_fft_input(x, N);

		auto* in = reinterpret_cast<T*>(&x_in[0]);

		std::vector<std::complex<T>> X(N);

		for (unsigned i = 0; i < N / 2; ++i)
		{
			X[i].real(in[2 * i]);
			X[i].imag(in[2 * i + 1]);
		}

		X = fft_(X, N / 2, NormalizationMode::backward);
		X.resize(N, { 0.0, 0.0 });

		auto* out = reinterpret_cast<T*>(&X[0]);

		unsigned nm1 = N - 1;
		unsigned nd2 = N / 2;
		unsigned n4 = (N / 4) - 1;

		for (unsigned i = 1; i <= n4; ++i)
		{
			auto im = nd2 - i;
			auto ip2 = i + nd2;
			auto ipm = im + nd2;
			out[2 * ip2] = (out[2 * i + 1] + out[2 * im + 1]) / 2;
			out[2 * ipm] = out[2 * ip2];
			out[2 * ip2 + 1] = -(out[2 * i] - out[2 * im]) / 2;
			out[2 * ipm + 1] = -out[2 * ip2 + 1];

			out[2 * i] = (out[2 * i] + out[2 * im]) / 2;
			out[2 * im] = out[2 * i];
			out[2 * i + 1] = (out[2 * i + 1] - out[2 * im + 1]) / 2;
			out[2 * im + 1] = -out[2 * i + 1];
		}

		out[2 * (N * 3) / 4] = out[2 * (N / 4) + 1];
		out[2 * nd2] = out[2 * 0 + 1];
		out[2 * (N * 3 / 4) + 1] = 0.0;
		out[2 * nd2 + 1] = 0.0;
		out[2 * (N / 4) + 1] = 0.0;
		out[2 * 0 + 1] = 0.0;

		auto l = static_cast<unsigned>(std::log2(N));
		unsigned le = 1 << l;
		unsigned le2 = le / 2;

		T u[2] = { 1.0, 0.0 };
		T s[2] = { static_cast<T>(std::cos(pi / le2)), static_cast<T>(-std::sin(pi / le2)) };
		T t[2];

		for (unsigned j = 1; j <= le2; ++j)
		{
			auto jm1 = j - 1;
			for (auto i = jm1; i <= nm1; i += le)
			{
				auto ip = i + le2;
				t[0] = out[2 * ip] * u[0] - out[2 * ip + 1] * u[1];
				t[1] = out[2 * ip] * u[1] + out[2 * ip + 1] * u[0];
				out[2 * ip] = out[2 * i] - t[0];
				out[2 * ip + 1] = out[2 * i + 1] - t[1];
				out[2 * i] += t[0];
				out[2 * i + 1] += t[1];
			}

			t[0] = u[0];
			u[0] = t[0] * s[0] - u[1] * s[1];
			u[1] = t[0] * s[1] + u[1] * s[0];
		}


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<T>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<T>(N); });
			break;
		}

		return X;
	}

	template<class T>
	std::vector<T> irfft_(const std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode)
	{
		if (x.empty()) return {};

		unsigned N = get_fft_length(x, n);

		auto x_in = resize_fft_input(x, N);

		for (unsigned k = N / 2 + 1; k < N; ++k)
		{
			x_in[k] = std::conj(x_in[N - k]);
		}
		for (auto& z : x_in)
		{
			z.real(z.real() + z.imag());
		}

		// Switch the mode because we are using the forward transform in the backward transform
		if (mode == NormalizationMode::backward)
		{
			mode = NormalizationMode::forward;
		}
		else if (mode == NormalizationMode::forward)
		{
			mode = NormalizationMode::backward;
		}
		auto xre = dsp::real(x_in);
		auto X = rfft_(xre, N, mode);

		for (unsigned i = 0; i < N; ++i)
		{
			X[i].real(X[i].real() + X[i].imag());
		}

		return dsp::real(X);
	}

#ifndef ZERO_DEPENDENCIES
	/* Wrapper functions for FFTW library functions of various precisions (high-performance for long signals) */
	auto fftw(const std::vector<std::complex<float>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<float>> X(N);
		std::vector<std::complex<float>> x_copy = x;
		x_copy.resize(N, 0.0);
		fftwf_complex* in = reinterpret_cast<fftwf_complex*>(&x_copy[0]);

		auto* out = reinterpret_cast<fftwf_complex*>(&X[0]);
		const auto p = fftwf_plan_dft_1d(N, in, out, sign, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<float>(N); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<float>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<float>(N); });
			}
			break;
		}

		fftwf_destroy_plan(p);
		return X;
	}
	auto fftw(const std::vector<std::complex<double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<double>> X(N);
		std::vector<std::complex<double>> x_copy = x;
		x_copy.resize(N, 0.0);
		fftw_complex* in = reinterpret_cast<fftw_complex*>(&x_copy[0]);

		auto* out = reinterpret_cast<fftw_complex*>(&X[0]);
		const auto p = fftw_plan_dft_1d(N, in, out, sign, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<double>(N); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<double>(N); });
			}
			break;
		}

		fftw_destroy_plan(p);
		return X;
	}
	auto fftw(const std::vector<std::complex<long double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<long double>> X(N);
		std::vector<std::complex<long double>> x_copy = x;
		x_copy.resize(N, 0.0);
		fftwl_complex* in = reinterpret_cast<fftwl_complex*>(&x_copy[0]);


		auto* out = reinterpret_cast<fftwl_complex*>(&X[0]);
		const auto p = fftwl_plan_dft_1d(N, in, out, sign, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<long double>(N); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<long double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [N](auto X) {return X / static_cast<long double>(N); });
			}
			break;
		}

		fftwl_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<double>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<double>> X(static_cast<size_t>(N / 2 + 1));
		std::vector<double> x_copy = x;
		x_copy.resize(N, 0.0);
		double* in = &x_copy[0];


		auto* out = reinterpret_cast<fftw_complex*>(&X[0]);
		auto* p = fftw_plan_dft_r2c_1d(N, in, out, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<double>(N); });
			break;
		}
		auto Xconj = dsp::conj(X);
		X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);
		fftw_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<float>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<float>> X(N);
		std::vector<float> x_copy = x;
		x_copy.resize(N, 0.0);
		float* in = &x_copy[0];


		auto* out = reinterpret_cast<fftwf_complex*>(&X[0]);
		const auto p = fftwf_plan_dft_r2c_1d(N, in, out, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<float>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<float>(N); });
			break;
		}

		auto Xconj = dsp::conj(X);
		X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);

		fftwf_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<long double>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<long double>> X(N);
		std::vector<long double> x_copy = x;
		x_copy.resize(N, 0.0);
		long double* in = &x_copy[0];

		auto* out = reinterpret_cast<fftwl_complex*>(&X[0]);
		const auto p = fftwl_plan_dft_r2c_1d(N, in, out, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<long double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [N](auto X) {return X / static_cast<long double>(N); });
			break;
		}

		auto Xconj = dsp::conj(X);
		X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);

		fftwl_destroy_plan(p);
		return X;
	}
	auto irfftw(const std::vector<std::complex<float>>& X, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (X.empty()) return std::vector<float>();

		unsigned N = get_fft_length(X, n);

		std::vector<float> x(N);
		std::vector<std::complex<float>> X_copy = X;
		X_copy.resize(N, 0.0);
		fftwf_complex* in = reinterpret_cast<fftwf_complex*>(&X_copy[0]);

		auto* out = &x[0];
		const auto p = fftwf_plan_dft_c2r_1d(N, in, out, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<float>(N); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<float>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			break;
		}

		fftwf_destroy_plan(p);
		return x;
	}
	auto irfftw(const std::vector<std::complex<double>>& X, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (X.empty()) return std::vector<double>();

		unsigned N = get_fft_length(X, n);

		std::vector<double> x(N);
		std::vector<std::complex<double>> X_copy = X;
		X_copy.resize(N, 0.0);
		fftw_complex* in = reinterpret_cast<fftw_complex*>(&X_copy[0]);


		auto* out = &x[0];
		const auto p = fftw_plan_dft_c2r_1d(N, in, out, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<double>(N); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			break;
		}

		fftw_destroy_plan(p);
		return x;
	}
	auto irfftw(const std::vector<std::complex<long double>>& X, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (X.empty()) return std::vector<long double>();

		unsigned N = get_fft_length(X, n);

		std::vector<long double> x(N);
		std::vector<std::complex<long double>> X_copy = X;
		X_copy.resize(N, 0.0);
		fftwl_complex* in = reinterpret_cast<fftwl_complex*>(&X_copy[0]);


		auto* out = &x[0];
		const auto p = fftwl_plan_dft_c2r_1d(N, in, out, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<long double>(N); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [N](auto x) {return x / static_cast<long double>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			break;
		}

		fftwl_destroy_plan(p);
		return x;
	}

#endif
	/// @endcond
} // .namespace dsp::fft



template<class T>
std::vector<std::complex<T>> dsp::fft::cfft(const std::vector<std::complex<T>>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend)
{
	switch (backend)
	{

	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (n > 100000)
		{
			return fftw(x, n, FFTW_FORWARD, FFTW_ESTIMATE, mode);
		}
#endif
		return fft_(x, n, mode);
	case backend::simple:
		return fft_(x, n, mode);
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return fftw(x, n, FFTW_FORWARD, FFTW_ESTIMATE, mode);
#else
		throw std::runtime_error("Library built without FFTW support!");
#endif
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<std::complex<T>> dsp::fft::icfft(const std::vector<std::complex<T>>& X, unsigned n, NormalizationMode mode, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (n > 100000)
		{
			return fftw(X, n, FFTW_BACKWARD, FFTW_ESTIMATE, mode);
		}
#endif
		return ifft_(X, n, mode);
	case backend::simple:
		return ifft_(X, n, mode);
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return fftw(X, n, FFTW_BACKWARD, FFTW_ESTIMATE, mode);
#else
		throw std::runtime_error("Library built without FFTW support!");
#endif
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<std::complex<T>> dsp::fft::rfft(const std::vector<T>& x, unsigned n, NormalizationMode mode, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (n > 100000)
		{
			return rfftw(x, n, FFTW_ESTIMATE, mode);
		}
#endif
		return rfft_(x, n, mode);
	case backend::simple:
		return rfft_(x, n, mode);
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return rfftw(x, n, FFTW_ESTIMATE, mode);
#else
		throw std::runtime_error("Library built without FFTW support!");
#endif
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<T> dsp::fft::irfft(const std::vector<std::complex<T>>& X, unsigned n, NormalizationMode mode, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (n > 100000)
		{
			return irfftw(X, n, FFTW_ESTIMATE, mode);
		}
#endif
		return irfft_(X, n, mode);
	case backend::simple:
		return irfft_(X, n, mode);
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return irfftw(X, n, FFTW_ESTIMATE, mode);
#else
		throw std::runtime_error("Library built without FFTW support!");
#endif
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<T> dsp::fft::logSquaredMagnitudeSpectrum(const std::vector<T>& signal, int N_fft,
	double relativeCutoff)
{
	auto spectrum = rfft(signal, 2 << (nextpow2(N_fft) - 1));

	const int finalFrequencyBinIdx = static_cast<const int>(relativeCutoff * static_cast<double>(spectrum.size()));

	auto logSquaredSpectrum = std::vector<T>(finalFrequencyBinIdx);


	std::transform(std::execution::par_unseq, spectrum.begin(), spectrum.begin() + finalFrequencyBinIdx, logSquaredSpectrum.begin(), &logSquaredMagnitude<T>);

	return logSquaredSpectrum;
}

template <class T>
std::vector<T> dsp::fft::fftconvolution(const std::vector<T>& volume, const std::vector<T>& kernel, convolution_mode mode)
{
	size_t size = static_cast<size_t>(std::pow(2, nextpow2(static_cast<unsigned>(volume.size() + kernel.size()) - 1)));
	std::vector<T> volume_padded(volume);
	volume_padded.resize(size);
	std::vector<T> kernel_padded(kernel);
	kernel_padded.resize(size);

	auto X = fft(volume_padded, static_cast<unsigned>(size));
	auto Y = fft(kernel_padded, static_cast<unsigned>(size));
	std::transform(X.begin(), X.end(), Y.begin(), X.begin(), std::multiplies<>());
	std::vector<T> result = ifft(X, static_cast<unsigned>(size));

	auto fullSize = volume.size() + kernel.size() - 1;

	switch (mode)
	{
	case convolution_mode::full:
		return { result.begin(), result.begin() + fullSize };
	case convolution_mode::valid:
		return centered<T>({ result.begin(), result.begin() + fullSize }, volume.size() - kernel.size() + 1);
	case convolution_mode::same:
		return centered<T>({ result.begin(), result.begin() + fullSize }, volume.size());
	default:
		throw std::runtime_error("Unknown convolution mode!");
		;
	}
}

template <class T>
std::vector<std::vector<T>> dsp::fft::spectrogram(const std::vector<T>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType)
{
	std::vector<std::vector<T>> spectrogram;

	// Algorithm to obtain the spectrogram is:
	// 0. Pre-emphasize the signal
	// 1. Create a number of overlapping frames from the signal
	// 2. Window each frame
	// 3. Calculate the log-squared-spectrum of each frame
	// -> Final result is a matrix of real values, each column representing one frame's log-squared magnitude spectrum

	// Pre-emphasis
	std::vector<T> b{ 1.0, static_cast<T>(-0.95) };
	std::vector<T> a{ 1.0 };
	auto preemph_signal = filter::filter(b, a, signal);

	// Split into frames
	auto frames = signalToFrames(signal, frameLength, static_cast<unsigned>(overlap_pct * frameLength));
	spectrogram.resize(frames.size());

	// Window each frame
	auto window = window::get_window<T>(windowType, frameLength);
	std::transform(std::execution::par_unseq,
		frames.begin(), frames.end(),
		frames.begin(),
		[window](auto& frame)
		{
			std::transform(frame.begin(), frame.end(), window.begin(), frame.begin(), std::multiplies<>());
			return frame;
		});

	// Calculate squared magnitude spectrum in dB
	const auto nFft = 2 << (nextpow2(frameLength) - 1);
	std::transform(std::execution::par_unseq,
		frames.begin(), frames.end(),
		spectrogram.begin(),
		[=](auto frame) {return logSquaredMagnitudeSpectrum<T>(frame, nFft, relativeCutoff); });

	return spectrogram;

}


template<class T>
std::vector<std::vector<T>> calcCosineBasisVectors(const unsigned int nBasisVectors)
{
	std::vector<std::vector<T>> basisVectors;
	basisVectors.reserve(nBasisVectors);

	T scalingFactorBn = static_cast<T>(std::sqrt(2.0 / nBasisVectors));
	T scalingFactorB1 = static_cast<T>(std::sqrt(1.0 / nBasisVectors));

	if (typeid(T) == typeid(float))
	{
		scalingFactorBn = std::sqrtf(static_cast<float>(2) / nBasisVectors);
		scalingFactorB1 = std::sqrtf(static_cast<float>(1) / nBasisVectors);
	}

	std::vector<T> b1(nBasisVectors, scalingFactorB1);
	basisVectors.push_back(b1);

	// i is the index of each basis vector.
	// n is the index of each element within each basis vector.
	for (unsigned int i = 1; i < nBasisVectors; ++i)
	{
		std::vector<T> bn;
		bn.reserve(nBasisVectors);

		for (unsigned int n = 0; n < nBasisVectors; ++n)
		{
			const T elem = static_cast<T>(cos((n + 0.5) * dsp::pi * i / nBasisVectors));  
			bn.push_back(elem * scalingFactorBn);
		}

		basisVectors.push_back(bn);
	}

	return basisVectors;
}


template<class T>
std::vector<T> dsp::fft::dct(std::vector<T>& signal, const unsigned int n, const dctType type)
{

	switch (type)
	{
	case dctType::dct1:
		throw std::runtime_error("DCT1 option is not implemented yet.");
	case dctType::dct3:
		throw std::runtime_error("DCT3 option is not implemented yet.");
	case dctType::dct4:
		throw std::runtime_error("DCT4 option is not implemented yet.");
	}
	// TODO: use static_assert() to catch not-implemented error during compile time.
	
	// Resize the signal. Zero-pads for signal lengths > n and
	// truncates for signal lengths < n.
	signal.resize(n);

	// Allocate output y.
	std::vector<T> y(n, 0);

	// Calculate B.
	std::vector<std::vector<T>> basisVectors = calcCosineBasisVectors<T>(n);

	// Calculate y.
	for (unsigned i = 0; i < n; ++i)
	{
		// Calculate inner product between the input signal and each basis vector (basis vectors are the rows in B^T).
		y[i] = std::inner_product(signal.begin(), signal.end(), basisVectors[i].begin(), static_cast<T>(0));

	}

	return y;
}


// Explicit template instantiation
template std::vector<std::complex<float>> dsp::fft::cfft(const std::vector<std::complex<float>>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<double>> dsp::fft::cfft(const std::vector<std::complex<double>>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<long double>> dsp::fft::cfft(const std::vector<std::complex<long double>>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);

template std::vector<std::complex<float>> dsp::fft::icfft(const std::vector<std::complex<float>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<double>> dsp::fft::icfft(const std::vector<std::complex<double>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<long double>> dsp::fft::icfft(const std::vector<std::complex<long double>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);

template std::vector<std::complex<float>> dsp::fft::rfft(const std::vector<float>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<double>> dsp::fft::rfft(const std::vector<double>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<std::complex<long double>> dsp::fft::rfft(const std::vector<long double>& x, unsigned n, dsp::fft::NormalizationMode mode, backend backend);

template std::vector<float> dsp::fft::irfft(const std::vector<std::complex<float>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<double> dsp::fft::irfft(const std::vector<std::complex<double>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);
template std::vector<long double> dsp::fft::irfft(const std::vector<std::complex<long double>>& X, unsigned n, dsp::fft::NormalizationMode mode, backend backend);

template std::vector<float> dsp::fft::fftconvolution(const std::vector<float>& volume, const std::vector<float>& kernel, convolution_mode mode);
template std::vector<double> dsp::fft::fftconvolution(const std::vector<double>& volume, const std::vector<double>& kernel, convolution_mode mode);
template std::vector<long double> dsp::fft::fftconvolution(const std::vector<long double>& volume, const std::vector<long double>& kernel, convolution_mode mode);

template std::vector<std::vector<float>> dsp::fft::spectrogram(const std::vector<float>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType);
template std::vector<std::vector<double>> dsp::fft::spectrogram(const std::vector<double>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType);
template std::vector<std::vector<long double>> dsp::fft::spectrogram(const std::vector<long double>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType);

template std::vector<float> dsp::fft::dct(std::vector<float>& signal, const unsigned int n, const dsp::fft::dctType type);
template std::vector<double> dsp::fft::dct(std::vector<double>& signal, const unsigned int n, const dsp::fft::dctType type);
template std::vector<long double> dsp::fft::dct(std::vector<long double>& signal, const unsigned int n, const dsp::fft::dctType type);

template std::vector<std::vector<float>> calcCosineBasisVectors(const unsigned int nBasisVectors);
template std::vector<std::vector<double>> calcCosineBasisVectors(const unsigned int nBasisVectors);
template std::vector<std::vector<long double>> calcCosineBasisVectors(const unsigned int nBasisVectors);
