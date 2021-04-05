#include "filter.h"

#include "dsp.h"

template <class T>
std::vector<T> dsp::filter::filter(std::vector<T> b,
	std::vector<T> a, const std::vector<T>& x)
{
	// Normalize filter coefficients
	auto a0 = a[0];
	for (auto& bi : b)
	{
		bi /= a0;
	}
	for (auto& ai : a)
	{
		ai /= a0;
	}

	if (a.size() == 1)
	{
		// This is the simple case without a denominator in the transfer function
		return dsp::convolve(x, b);
	}

	// This is the more complicated case with both numerator and denominator

	std::vector<T> y(x.size());  // The filtered output
	auto n = std::max(b.size(), a.size());
	a.resize(n);
	b.resize(n);
	std::vector<T> w_m(n + 1, 0);  // w(m)
	std::vector<T> w_mm1(n + 1, 0);  // w(m-1)
	for (size_t m = 0; m < x.size(); ++m)
	{
		y[m] = b[0] * x[m] + w_mm1[1];
		for (size_t k = n-1; k > 0; --k)
		{
			w_m[k] = b[k] * x[m] + w_mm1[k + 1] /* zero for k = n-1 */ - a[k] * y[m];
		}		
		w_mm1 = w_m;
	}

	return y;
}

template<class T>
std::vector<T> dsp::filter::lpc(const std::vector<T>& x, unsigned N)
{
	std::vector<T> r;
	auto numCoefficients = N + 1;
	r.resize(numCoefficients);

	for (size_t i = 0; i <= N; i++)
	{
		r[i] = 0.0;
		for (size_t j = 0; j < x.size() - i; ++j)
		{
			r[i] += x[j] * x[j + i];
		}
	}

	// Levinson-Durbin
	T E = r[0];
	std::vector<T> alpha(numCoefficients, 0);
	alpha[0] = 1.0;
	std::vector<T> beta(numCoefficients, 0);
	std::vector<T> z(numCoefficients, 0);

	for (unsigned p = 1; p <= N; ++p)
	{
		T q = 0.0;
		for (unsigned i = 0; i < p; ++i) { q += alpha[i] * r[p - i]; }
		if (E == 0.0) { E = static_cast<T>(0.0001); }
		z[p] = -q / E;
		alpha[p] = 0.0;
		for (unsigned i = 0; i <= p; ++i) { beta[i] = alpha[i] + z[p] * alpha[p - i]; }
		for (unsigned i = 0; i <= p; ++i) { alpha[i] = beta[i]; }

		E = E * (static_cast<T>(1.0) - z[p] * z[p]);
	}

	std::vector<T> coeff;
	coeff.resize(N + 1);
	coeff[0] = 1;
	for (unsigned i = 1; i <= N; ++i) { coeff[i] = -alpha[i]; }

	return coeff;
}

template <class T>
std::vector<T> dsp::filter::medianfilter(const std::vector<T>& x, size_t kernel_size)
{
	if (kernel_size % 2 != 1) { throw std::runtime_error("Kernel size should be odd!"); }

	// Pad vector with zeros at the edges
	std::vector<T> padded_x = pad(x, { kernel_size / 2, kernel_size / 2 });
	std::vector<T> out;

	// Find median in each kernel
	for (size_t k = 0; k < x.size(); ++k)
	{
		out.push_back(median<T>(padded_x.begin() + k, padded_x.begin() + k + kernel_size));
	}
	
	return out;
}

template <class T>
dsp::Signal<T> dsp::filter::medianfilter(const Signal<T>& x, typename Signal<T>::size_type kernel_size)
{
	return Signal(x.getSamplingRate_Hz(), medianfilter<T>(x.getSamples(), kernel_size));
}


// Explicit template instantiation
template std::vector<float> dsp::filter::filter(std::vector<float> b,
	std::vector<float> a, const std::vector<float>& x);
template std::vector<double> dsp::filter::filter(std::vector<double> b,
	std::vector<double> a, const std::vector<double>& x);
template std::vector<long double> dsp::filter::filter(std::vector<long double> b,
	std::vector<long double> a, const std::vector<long double>& x);

template std::vector<float> dsp::filter::lpc(const std::vector<float>& x, unsigned N);
template std::vector<double> dsp::filter::lpc(const std::vector<double>& x, unsigned N);
template std::vector<long double> dsp::filter::lpc(const std::vector<long double>& x, unsigned N);

template std::vector<float> dsp::filter::medianfilter(const std::vector<float>& x, size_t kernel_size);
template std::vector<double> dsp::filter::medianfilter(const std::vector<double>& x, size_t kernel_size);
template std::vector<long double> dsp::filter::medianfilter(const std::vector<long double>& x, size_t kernel_size);

template dsp::Signal<float> dsp::filter::medianfilter(const dsp::Signal<float>& x, Signal<float>::size_type kernel_size);
template dsp::Signal<double> dsp::filter::medianfilter(const dsp::Signal<double>& x, Signal<double>::size_type kernel_size);
template dsp::Signal<long double> dsp::filter::medianfilter(const dsp::Signal<long double>& x, Signal<long double>::size_type kernel_size);