#include "fft.h"

#include <fftw3.h>

namespace dsp::fft
{
	/// <summary>
	/// Helper function to get the FFT length
	/// </summary>
	template<class T>
	auto get_fft_length(const std::vector<T>& x, unsigned n)
	{
		if (n == 0)
		{
			return static_cast<unsigned>(x.size());
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
	std::vector<std::complex<T>> fft_(std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return {};

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<T>> X = resize_fft_input(x, N);
		auto* in = reinterpret_cast<T*>(&X[0]);

		int nm1 = N - 1;
		int nd2 = N / 2;
		unsigned j = nd2;
		for (unsigned i = 1; i <= N - 2; ++i)
		{
			if (i < j)
			{
				T t[2] = { in[2 * j], in[2 * j + 1] };
				in[2 * j] = in[2 * i];
				in[2 * j + 1] = in[2 * i + 1];
				in[2 * i] = t[0];
				in[2 * i + 1] = t[1];
			}
			unsigned k = nd2;

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
	std::vector<std::complex<T>> ifft_(std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
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
		X = fft_(X, N, mode, overwrite_x);

		return dsp::conj(X);
	}

	template<class T>
	std::vector<std::complex<T>> rfft_(std::vector<T>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
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

		X = fft_(X, N / 2, NormalizationMode::backward, overwrite_x);
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
	std::vector<T> irfft_(std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
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
		auto X = rfft_(xre, N, mode, overwrite_x);

		for (unsigned i = 0; i < N; ++i)
		{
			X[i].real(X[i].real() + X[i].imag());
		}

		return dsp::real(X);
	}


	/* Wrapper functions for FFTW library functions of various precisions (high-performance for long signals) */
	auto fftw(std::vector<std::complex<float>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<float>> X(N);
		std::vector<std::complex<float>> x_copy;
		fftwf_complex* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = reinterpret_cast<fftwf_complex*>(&x[0]);
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = reinterpret_cast<fftwf_complex*>(&x_copy[0]);
		}

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
	auto fftw(std::vector<std::complex<double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<double>> X(N);
		std::vector<std::complex<double>> x_copy;
		fftw_complex* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = reinterpret_cast<fftw_complex*>(&x[0]);
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = reinterpret_cast<fftw_complex*>(&x_copy[0]);
		}

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
	auto fftw(std::vector<std::complex<long double>>& x, unsigned n, int sign, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<long double>> X(N);
		std::vector<std::complex<long double>> x_copy;
		fftwl_complex* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = reinterpret_cast<fftwl_complex*>(&x[0]);
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = reinterpret_cast<fftwl_complex*>(&x_copy[0]);
		}

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
	auto rfftw(std::vector<double>& x, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<double>> X(static_cast<size_t>(N / 2 + 1));
		std::vector<double> x_copy;
		double* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = &x[0];
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = &x_copy[0];
		}

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
	auto rfftw(std::vector<float>& x, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<float>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<float>> X(N);
		std::vector<float> x_copy;
		float* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = &x[0];
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = &x_copy[0];
		}

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
	auto rfftw(std::vector<long double>& x, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<long double>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<long double>> X(N);
		std::vector<long double> x_copy;
		long double* in;
		if (overwrite_x)
		{
			x.resize(N, 0.0);
			in = &x[0];
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, 0.0);
			in = &x_copy[0];
		}

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
	auto irfftw(std::vector<std::complex<float>>& X, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_X)
	{
		if (X.empty()) return std::vector<float>();

		unsigned N = get_fft_length(X, n);

		std::vector<float> x(N);
		std::vector<std::complex<float>> X_copy;
		fftwf_complex* in;
		if (overwrite_X)
		{
			X.resize(N, 0.0);
			in = reinterpret_cast<fftwf_complex*>(&X[0]);
		}
		else
		{
			X_copy = X;
			X_copy.resize(N, 0.0);
			in = reinterpret_cast<fftwf_complex*>(&X_copy[0]);
		}

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
	auto irfftw(std::vector<std::complex<double>>& X, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_X)
	{
		if (X.empty()) return std::vector<double>();

		unsigned N = get_fft_length(X, n);

		std::vector<double> x(N);
		std::vector<std::complex<double>> X_copy;
		fftw_complex* in;
		if (overwrite_X)
		{
			X.resize(N, 0.0);
			in = reinterpret_cast<fftw_complex*>(&X[0]);
		}
		else
		{
			X_copy = X;
			X_copy.resize(N, 0.0);
			in = reinterpret_cast<fftw_complex*>(&X_copy[0]);
		}

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
	auto irfftw(std::vector<std::complex<long double>>& X, unsigned n, unsigned flags, NormalizationMode mode, bool overwrite_X)
	{
		if (X.empty()) return std::vector<long double>();

		unsigned N = get_fft_length(X, n);

		std::vector<long double> x(N);
		std::vector<std::complex<long double>> X_copy;
		fftwl_complex* in;
		if (overwrite_X)
		{
			X.resize(N, 0.0);
			in = reinterpret_cast<fftwl_complex*>(&X[0]);
		}
		else
		{
			X_copy = X;
			X_copy.resize(N, 0.0);
			in = reinterpret_cast<fftwl_complex*>(&X_copy[0]);
		}

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
}

template<class T>
std::vector<std::complex<T>> dsp::fft::fft(std::vector<std::complex<T>>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend)
{
	switch (backend)
	{

	case backend::automatic:
		if (n > 100000)
		{
			return fftw(x, n, FFTW_FORWARD, FFTW_ESTIMATE, mode, overwrite_x);
		}
		return fft_(x, n, mode, overwrite_x);
	case backend::simple:
		return fft_(x, n, mode, overwrite_x);
	case backend::fftw:
		return fftw(x, n, FFTW_FORWARD, FFTW_ESTIMATE, mode, overwrite_x);
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<std::complex<T>> dsp::fft::ifft(std::vector<std::complex<T>>& X, unsigned n, NormalizationMode mode, bool overwrite_X, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
		if (n > 100000)
		{
			return fftw(X, n, FFTW_BACKWARD, FFTW_ESTIMATE, mode, overwrite_X);
		}
		return ifft_(X, n, mode, overwrite_X);
	case backend::simple:
		return ifft_(X, n, mode, overwrite_X);
	case backend::fftw:
		return fftw(X, n, FFTW_BACKWARD, FFTW_ESTIMATE, mode, overwrite_X);
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<std::complex<T>> dsp::fft::rfft(std::vector<T>& x, unsigned n, NormalizationMode mode, bool overwrite_x, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
		if (n > 100000)
		{
			return rfftw(x, n, FFTW_ESTIMATE, mode, overwrite_x);
		}
		return rfft_(x, n, mode, overwrite_x);
	case backend::simple:
		return rfft_(x, n, mode, overwrite_x);
	case backend::fftw:
		return rfftw(x, n, FFTW_ESTIMATE, mode, overwrite_x);
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

template <class T>
std::vector<T> dsp::fft::irfft(std::vector<std::complex<T>>& X, unsigned n, NormalizationMode mode, bool overwrite_X, backend backend)
{
	switch (backend)
	{
	case backend::automatic:
		if (n > 100000)
		{
			return irfftw(X, n, FFTW_ESTIMATE, mode, overwrite_X);
		}
		return irfft_(X, n, mode, overwrite_X);
	case backend::simple:
		return irfft_(X, n, mode, overwrite_X);
	case backend::fftw:
		return irfftw(X, n, FFTW_ESTIMATE, mode, overwrite_X);
	default:
		throw std::runtime_error("Unknown backend selected!");
	}
}

// Explicit template instantiation
template std::vector<std::complex<float>> dsp::fft::fft(std::vector<std::complex<float>>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<double>> dsp::fft::fft(std::vector<std::complex<double>>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<long double>> dsp::fft::fft(std::vector<std::complex<long double>>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<float>> dsp::fft::ifft(std::vector<std::complex<float>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<double>> dsp::fft::ifft(std::vector<std::complex<double>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<long double>> dsp::fft::ifft(std::vector<std::complex<long double>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<float>> dsp::fft::rfft(std::vector<float>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<double>> dsp::fft::rfft(std::vector<double>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<std::complex<long double>> dsp::fft::rfft(std::vector<long double>& x, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_x, backend backend);
template std::vector<float> dsp::fft::irfft(std::vector<std::complex<float>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_X, backend backend);
template std::vector<double> dsp::fft::irfft(std::vector<std::complex<double>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_X, backend backend);
template std::vector<long double> dsp::fft::irfft(std::vector<std::complex<long double>>& X, unsigned n, dsp::fft::NormalizationMode mode, bool overwrite_X, backend backend);