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

	/* Straight-forward FFT implementations (high-performance for short signals) */
	template<class T>
	std::vector<std::complex<T>> fft_(std::vector<std::complex<T>>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return std::vector<std::complex<T>>();

		unsigned N = get_fft_length(x, n);

		std::vector<std::complex<T>> x_copy;
		std::complex<T>* in;
		
		if (overwrite_x)
		{
			x.resize(N, { 0.0, 0.0 });
			in = &x[0];
		}
		else
		{
			x_copy = x;
			x_copy.resize(N, { 0.0, 0.0 });
			in = &x_copy[0];
		}
		int nm1 = N - 1;
		int nd2 = N / 2;
		unsigned j = nd2;
		for(unsigned i = 1; i <= N-2; ++i)
		{
			if(i < j)
			{
				std::swap(in[j], in[i]);
			}
			unsigned k = nd2;

			while(k <= j)
			{
				j -= k;
				k /= 2;
			}
			j += k;
		}
		auto exponent = static_cast<unsigned>(std::log2(n));
		for (unsigned l = 1; l <= exponent; ++l)
		{
			auto le = 1 << l;
			auto le2 = le / 2;
			std::complex<T> u(1.0, 0.0);
			std::complex<T> s(static_cast<T>(std::cos(pi / le2)), static_cast<T>(-1 * std::sin(pi / le2)));

			for (auto j = 1; j <= le2; ++j)
			{
				auto jm1 = j - 1;
				std::complex<T> t;
				for (auto i = jm1; i <= nm1; i += le)
				{
					auto ip = i + le2;
					t = std::complex<T>(in[ip].real() * u.real() - in[ip].imag() * u.imag(), in[ip].real() * u.imag() + in[ip].imag() * u.real());
					in[ip] = in[i] - t;
					in[i] += t;
				}

				t.real(u.real());
				u = std::complex<T>(t.real() * s.real() - u.imag() * s.imag(), t.real() * s.imag() + u.imag() * s.real());			
			}
		}

		std::vector<std::complex<T>>* X;
		if (overwrite_x)
		{
			X = &x;
		}
		else
		{
			X = &x_copy;
		}
		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(X->begin(), X->end(),
				X->begin(), [N](auto X) {return X / static_cast<T>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(X->begin(), X->end(),
				X->begin(), [N](auto X) {return X / static_cast<T>(N); });
			break;
		}

		return *X;
	}

	// TODO: ifft_

	template<class T>
	auto resize_fft_input(std::vector<T> x, unsigned n)
	{
		x.resize(n, T());
		return x;
	}
	
	template<class T>
	std::vector<std::complex<T>> rfft_(std::vector<T>& x, unsigned n, NormalizationMode mode, bool overwrite_x)
	{
		if (x.empty()) return {};

		unsigned N = get_fft_length(x, n);

		std::vector<T> in = resize_fft_input(x, n);

		std::vector<std::complex<T>> out(N);

		for(unsigned i = 0; i < N/2; ++i)
		{
			out[i].real(in[2 * i]);
			out[i].imag(in[2 * i + 1]);
		}

		out = fft_(out, N / 2, NormalizationMode::backward, true);
		out.resize(N, { 0.0, 0.0 });

		unsigned nm1 = N - 1;
		unsigned nd2 = N / 2;
		unsigned n4 = (N / 4) - 1;

		for (unsigned i = 1; i <= n4; ++i)
		{
			auto im = nd2 - i;
			auto ip2 = i + nd2;
			auto ipm = im + nd2;
			out[ip2].real((out[i].imag() + out[im].imag()) / 2);
			out[ipm].real(out[ip2].real());
			out[ip2].imag(-(out[i].real() - out[im].real()) / 2);
			out[ipm].imag(-out[ip2].imag());

			out[i].real((out[i].real() + out[im].real()) / 2);
			out[im].real(out[i].real());
			out[i].imag((out[i].imag() - out[im].imag()) / 2);
			out[im].imag(-out[i].imag());
		}

		out[(N * 3) / 4].real(out[N / 4].imag());
		out[nd2].real(out[0].imag());
		out[(N * 3 / 4)].imag(0.0);
		out[nd2].imag(0.0);
		out[N / 4].imag(0.0);
		out[0].imag(0.0);

		auto l = static_cast<unsigned>(std::log2(N));
		unsigned le = 1 << l;
		unsigned le2 = le / 2;

		std::complex<T> u(1.0, 0.0);
		std::complex<T> s(static_cast<T>(std::cos(pi / le2)), static_cast<T>(-std::sin(pi / le2)));
		std::complex<T> t;

		for (unsigned j = 1; j <= le2; ++j)
		{
			auto jm1 = j - 1;
			for (auto i = jm1; i <= nm1; i+= le)
			{
				auto ip = i + le2;
				t = { out[ip].real() * u.real() - out[ip].imag() * u.imag(), out[ip].real() * u.imag() + out[ip].imag() * u.real() };
				out[ip] = out[i] - t;
				out[i] += t;
			}

			t.real(u.real());
			u.real(t.real() * s.real() - u.imag() * s.imag());
			u.imag(t.real() * s.imag() + u.imag() * s.real());
		}
		
		
		switch (mode)
		{
		case NormalizationMode::backward:
			break;
		case NormalizationMode::ortho:
			std::transform(out.begin(), out.end(),
				out.begin(), [N](auto X) {return X / static_cast<T>(sqrt(N)); });
			break;
		case NormalizationMode::forward:
			std::transform(out.begin(), out.end(),
				out.begin(), [N](auto X) {return X / static_cast<T>(N); });
			break;
		}

		return out;
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

		std::vector<std::complex<double>> X(static_cast<size_t>(N/2+1));
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
	switch(backend)
	{
		
	case backend::automatic:
		if(n > 10000000)
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
	return fftw(X, n, FFTW_BACKWARD, FFTW_ESTIMATE, mode, overwrite_X);
}

template <class T>
std::vector<std::complex<T>> dsp::fft::rfft(std::vector<T>& x, unsigned n, NormalizationMode mode, bool overwrite_x, backend backend)
{
	switch(backend)
	{
	case backend::automatic:
		if(n > 1000000)
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
	return irfftw(X, n, FFTW_ESTIMATE, mode, overwrite_X);
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