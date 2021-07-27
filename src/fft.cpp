#include "fft.h"

#include <algorithm>
#include <execution>

#ifndef ZERO_DEPENDENCIES
#include <fftw3.h>
#endif

#include "filter.h"
#include "Signal.h"


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
			return static_cast<unsigned> (2 << (nextpow2(n) - 1));
		}
		
		return n;
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

		// DEBUG
		auto N = n;

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

	template <class T>
	std::vector<std::complex<T>> dft_(const std::vector<T>& x, unsigned n, NormalizationMode mode)
	{
		std::vector<std::complex<T>> X;		

		for (unsigned k = 0; k <= n/2; ++k)
		{
			T re{ 0.0 };
			T im{ 0.0 };

			for (unsigned i = 0; i < n; ++i)
			{
				double angle = (2.0 * pi * k * i) / static_cast<double>(n);
				re += static_cast<T>(x[i] * std::cos(angle));
				im -= static_cast<T>(x[i] * std::sin(angle)); 
			}

			if (mode == NormalizationMode::forward)
			{
				im = -im;
				re /= static_cast<T>(n / 2);
				im /= static_cast<T>(n / 2);
				if ((k == 0) || (k == n/2)) { re /= 2; }
			}
			X.emplace_back(re, im);
		}
		auto Xconj = dsp::conj(X);
		if (n % 2 == 0)
		{
			X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);
		}
		else
		{
			X.insert(X.end(), Xconj.rbegin(), Xconj.rend() - 1);
		}
		return X;
	}

	template <class T>
	std::vector<T> idft_(const std::vector<std::complex<T>> X, unsigned n, NormalizationMode mode)
	{
		T re{ 0.0 }, im{ 0.0 };
		
		std::vector<T> x(n, 0.0);
		
		for (unsigned k = 0; k <= n / 2; ++k)
		{
			if (mode == NormalizationMode::backward)
			{
				im = -X[k].imag() / static_cast<T>(n / 2);
				re = X[k].real() / static_cast<T>(n / 2);
				if ((k == 0) || (k == n / 2)) { re /= 2.0; }
			}
			else
			{
				re = X[k].real();
				im = X[k].imag();
			}

			for (unsigned i = 0; i < n; ++i)
			{
				T angle = static_cast<T>(2.0 * pi * k * i) / static_cast<T>(n);
				x[i] += re * std::cos(angle) + im * std::sin(angle);
			}
		}
		return x;
	}

#ifndef ZERO_DEPENDENCIES
	/* Wrapper functions for FFTW library functions of various precisions (high-performance for long signals) */
	auto fftw(const std::vector<std::complex<float>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		std::vector<std::complex<float>> X(n);
		std::vector<std::complex<float>> x_copy = x;
		x_copy.resize(n, 0.0);
		fftwf_complex* in = reinterpret_cast<fftwf_complex*>(&x_copy[0]);

		auto* out = reinterpret_cast<fftwf_complex*>(&X[0]);
		const auto p = fftwf_plan_dft_1d(n, in, out, sign, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<float>(n); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<float>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<float>(n); });
			}
			break;
		}

		fftwf_destroy_plan(p);
		return X;
	}
	auto fftw(const std::vector<std::complex<double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		std::vector<std::complex<double>> X(n);
		std::vector<std::complex<double>> x_copy = x;
		x_copy.resize(n, 0.0);
		fftw_complex* in = reinterpret_cast<fftw_complex*>(&x_copy[0]);

		auto* out = reinterpret_cast<fftw_complex*>(&X[0]);
		const auto p = fftw_plan_dft_1d(n, in, out, sign, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<double>(n); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<double>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<double>(n); });
			}
			break;
		}

		fftw_destroy_plan(p);
		return X;
	}
	auto fftw(const std::vector<std::complex<long double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		std::vector<std::complex<long double>> X(n);
		std::vector<std::complex<long double>> x_copy = x;
		x_copy.resize(n, 0.0);
		fftwl_complex* in = reinterpret_cast<fftwl_complex*>(&x_copy[0]);


		auto* out = reinterpret_cast<fftwl_complex*>(&X[0]);
		const auto p = fftwl_plan_dft_1d(n, in, out, sign, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			if (sign == FFTW_BACKWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<long double>(n); });
			}
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<long double>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			if (sign == FFTW_FORWARD)
			{
				std::transform(X.begin(), X.end(),
					X.begin(), [n](auto X) {return X / static_cast<long double>(n); });
			}
			break;
		}

		auto Xconj = dsp::conj(X);
		if (n % 2 == 0)
		{
			X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);
		}
		else
		{
			X.insert(X.end(), Xconj.rbegin(), Xconj.rend() - 1);
		}
		fftwl_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<double>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		std::vector<std::complex<double>> X(static_cast<size_t>(n / 2 + 1));
		std::vector<double> x_copy = x;
		x_copy.resize(n, 0.0);
		double* in = &x_copy[0];


		auto* out = reinterpret_cast<fftw_complex*>(&X[0]);
		auto* p = fftw_plan_dft_r2c_1d(n, in, out, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<double>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<double>(n); });
			break;
		}
		auto Xconj = dsp::conj(X);
		if (n % 2 == 0)
		{
			X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);
		}
		else
		{
			X.insert(X.end(), Xconj.rbegin(), Xconj.rend() - 1);
		}
		fftw_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<float>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		std::vector<std::complex<float>> X(static_cast<size_t>(n / 2 + 1));
		std::vector<float> x_copy = x;
		x_copy.resize(n, 0.0);
		float* in = &x_copy[0];


		auto* out = reinterpret_cast<fftwf_complex*>(&X[0]);
		const auto p = fftwf_plan_dft_r2c_1d(n, in, out, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<float>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<float>(n); });
			break;
		}

		auto Xconj = dsp::conj(X);
		if (n % 2 == 0)
		{
			X.insert(X.end(), Xconj.rbegin() + 1, Xconj.rend() - 1);
		}
		else
		{
			X.insert(X.end(), Xconj.rbegin(), Xconj.rend() - 1);
		}
		fftwf_destroy_plan(p);
		return X;
	}
	auto rfftw(const std::vector<long double>& x, unsigned n, unsigned flags, NormalizationMode mode)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		std::vector<std::complex<long double>> X(static_cast<size_t>(n / 2 + 1));
		std::vector<long double> x_copy = x;
		x_copy.resize(n, 0.0);
		long double* in = &x_copy[0];

		auto* out = reinterpret_cast<fftwl_complex*>(&X[0]);
		const auto p = fftwl_plan_dft_r2c_1d(n, in, out, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<long double>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			std::transform(X.begin(), X.end(),
				X.begin(), [n](auto X) {return X / static_cast<long double>(n); });
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

		std::vector<float> x(n);
		std::vector<std::complex<float>> X_copy = X;
		X_copy.resize(n, 0.0);
		fftwf_complex* in = reinterpret_cast<fftwf_complex*>(&X_copy[0]);

		auto* out = &x[0];
		const auto p = fftwf_plan_dft_c2r_1d(n, in, out, flags);

		fftwf_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<float>(n); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<float>(sqrt(n)); });
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

		std::vector<double> x(n);
		std::vector<std::complex<double>> X_copy = X;
		X_copy.resize(n, 0.0);
		fftw_complex* in = reinterpret_cast<fftw_complex*>(&X_copy[0]);


		auto* out = &x[0];
		const auto p = fftw_plan_dft_c2r_1d(n, in, out, flags);

		fftw_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<double>(n); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<double>(sqrt(n)); });
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

		std::vector<long double> x(n);
		std::vector<std::complex<long double>> X_copy = X;
		X_copy.resize(n, 0.0);
		fftwl_complex* in = reinterpret_cast<fftwl_complex*>(&X_copy[0]);


		auto* out = &x[0];
		const auto p = fftwl_plan_dft_c2r_1d(n, in, out, flags);

		fftwl_execute(p);


		switch (mode)
		{
		case NormalizationMode::backward:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<long double>(n); });
			break;
		case NormalizationMode::ortho:
			std::transform(x.begin(), x.end(),
				x.begin(), [n](auto x) {return x / static_cast<long double>(sqrt(n)); });
			break;
		case NormalizationMode::forward:
			break;
		}

		fftwl_destroy_plan(p);
		return x;
	}

#endif
	/// @endcond
}



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
	auto N = get_fft_length(x, n);
	
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (N > 100000)
		{
			return rfftw(x, N, FFTW_ESTIMATE, mode);
		}
#endif
		if (ispow2(N))
		{
			return rfft_(x, N, mode);
		}
		else
		{
			return dft_(x, N, mode);
		}
	case backend::simple:
		if (ispow2(N))
		{
			return rfft_(x, N, mode);
		}
		else
		{			
			return dft_(x, N, mode);
		}
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return rfftw(x, N, FFTW_ESTIMATE, mode);
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
	auto N = get_fft_length(X, n);
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (N > 100000)
		{
			return irfftw(X, N, FFTW_ESTIMATE, mode);
		}
#endif
		if (ispow2(N))
		{
			return irfft_(X, N, mode);
		}
		else
		{
			return idft_(X, N, mode);
		}
	case backend::simple:
		if (ispow2(N))
		{
			return irfft_(X, N, mode);
		}
		else
		{
			return idft_(X, N, mode);
		}
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		return irfftw(X, N, FFTW_ESTIMATE, mode);
#else
		throw std::runtime_error("Library built without FFTW support!");
#endif
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template<class T>
std::vector<std::complex<T>> dsp::fft::dft(const std::vector<T>& x, unsigned n, NormalizationMode mode, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
#ifndef ZERO_DEPENDENCIES
		if (n > 100000)
		{
			//return dftw(X, n, FFTW_ESTIMATE, mode);
			throw std::runtime_error("Not implemented yet!");
		}
#endif
		return dft_(x, n, mode);
	case backend::simple:
		return dft_(x, n, mode);
	case backend::fftw:
#ifndef ZERO_DEPENDENCIES
		//return dftw(X, n, FFTW_ESTIMATE, mode);
		throw std::runtime_error("Not implemented yet!");
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
	auto spectrum = rfft(signal, 2 << (nextpow2(N_fft)-1));

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

template<class T>
std::vector<std::vector<std::complex<T>>> dsp::fft::stft(const std::vector<T>& x, unsigned frameLength, int overlap, window::type window, unsigned fftLength)
{
	// Algorithm:
	// 1. Split input into frames of frameLength samples, overlapping by overlap samples. Pad with zeros at the end.
	// 2. Window each frame using the specified window type.
	// 3. Calculate the FFT of every frame

	auto frames = signalToFrames(x, frameLength, overlap, false);

	// TODO: Figure out how to allow windows with additional parameters
	auto w = dsp::window::get_window<T>(window, frameLength, false);

	for (auto& frame : frames)
	{
		std::transform(frame.begin(), frame.end(), w.begin(), frame.begin(), std::multiplies<>());
	}


	std::vector<std::vector<std::complex<T>>> stft;
	stft.reserve(frames.size());
	for (const auto& frame : frames)
	{
		stft.emplace_back(dsp::fft::rfft(frame, fftLength));
	}
	
	return stft;
}


template <class T>
std::vector<std::vector<T>> dsp::fft::spectrogram(const std::vector<T>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType, bool logSquared)
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

	auto N = 2 << (nextpow2(frameLength) - 1);

	auto overlap = static_cast<int>(overlap_pct * frameLength);

	auto stft = dsp::fft::stft(preemph_signal, frameLength, overlap, windowType, N);

	const int finalFrequencyBinIdx = std::min(N, static_cast<const int>(relativeCutoff * static_cast<double>(N)) + 1);

	spectrogram.reserve(stft.size());

	for (const auto& frame : stft)
	{
		auto magnitudeSpectrum = std::vector<T>(finalFrequencyBinIdx);
		if (logSquared)
		{
			std::transform(std::execution::par_unseq, frame.begin(), frame.begin() + finalFrequencyBinIdx, magnitudeSpectrum.begin(), &logSquaredMagnitude<T>);
		}
		else
		{
			std::transform(std::execution::par_unseq, frame.begin(), frame.begin() + finalFrequencyBinIdx, magnitudeSpectrum.begin(), &std::abs<T>);
		}
		spectrogram.push_back(magnitudeSpectrum);
	}

	return spectrogram;

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

template std::vector<std::vector<std::complex<float>>> dsp::fft::stft(const std::vector<float> & x, unsigned frameLength, int overlap, window::type window, unsigned fftLength);
template std::vector<std::vector<std::complex<double>>> dsp::fft::stft(const std::vector<double> & x, unsigned frameLength, int overlap, window::type window, unsigned fftLength);
template std::vector<std::vector<std::complex<long double>>> dsp::fft::stft(const std::vector<long double> & x, unsigned frameLength, int overlap, window::type window, unsigned fftLength);

template std::vector<std::vector<float>> dsp::fft::spectrogram(const std::vector<float>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType, bool logSquared);
template std::vector<std::vector<double>> dsp::fft::spectrogram(const std::vector<double>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType, bool logSquared);
template std::vector<std::vector<long double>> dsp::fft::spectrogram(const std::vector<long double>& signal, unsigned frameLength,
	double overlap_pct, int samplingRate, double relativeCutoff, window::type windowType, bool logSquared);

template std::vector<float> dsp::fft::logSquaredMagnitudeSpectrum(const std::vector<float>& signal, int N_fft,
	double relativeCutoff);
template std::vector<double> dsp::fft::logSquaredMagnitudeSpectrum(const std::vector<double>& signal, int N_fft,
	double relativeCutoff);
template std::vector<long double> dsp::fft::logSquaredMagnitudeSpectrum(const std::vector<long double>& signal, int N_fft,
	double relativeCutoff);