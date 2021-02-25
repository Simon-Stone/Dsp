#include "Window.h"

#include <algorithm>

#include "Signal.h"

namespace Dsp::Window
{
	template<class T>
	std::vector<T> boxcar(unsigned N, bool sym)
	{
		if (N <= 0) return std::vector<T>();

		return std::vector<T>(N, static_cast<T>(1));
	}

	template<class T>
	std::vector<T> triang(unsigned N, bool sym)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);
		auto n = Dsp::arange<T>(1.0, static_cast<T>((M + 1) / 2 + 1));

		std::vector<T> w;
		w.resize(n.size(), 0.0);
		if (M % 2 == 0)
		{
			std::transform(n.begin(), n.end(),
				w.begin(),
				[M](auto n) {return static_cast<T>((2 * n - 1.0) / M); });
			w.insert(w.end(), w.rbegin(), w.rend());
		}
		else
		{
			std::transform(n.begin(), n.end(),
				w.begin(),
				[M](auto n) {return static_cast<T>(2 * n / (M + 1.0)); });
			w.insert(w.end(), w.rbegin() + 1, w.rend());
		}

		return Utilities::truncate(w, needs_trunc);
	}

	template<class T>
	std::vector<T> general_cosine(unsigned N, const std::vector<T>& a, bool sym)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);

		auto fac = Dsp::linspace<T>(static_cast<T>(-pi), static_cast<T>(pi), M);
		std::vector<T> w(M, 0);
		for (size_t k = 0; k < a.size(); ++k)
		{
			auto ak = a[k];
			std::transform(w.begin(), w.end(),
				fac.begin(),
				w.begin(),
				[k, ak](auto w, auto fac) {return w + ak * std::cos(k * fac); });
		}
		return Utilities::truncate(w, needs_trunc);
	}

	template <class T>
	std::vector<T> general_hamming(unsigned N, double alpha, bool sym)
	{
		return general_cosine<T>(N, { static_cast<T>(alpha), static_cast<T>(1.0 - alpha) }, sym);
	}

	template<class T>
	std::vector<T> blackman(unsigned N, bool sym)
	{
		return general_cosine<T>(N, { static_cast<T>(0.42), static_cast<T>(0.50), static_cast<T>(0.08) }, sym);
	}

	template <class T>
	std::vector<T> hamming(unsigned N, bool sym)
	{
		return general_hamming<T>(N, 0.54, sym);
	}

	template <class T>
	std::vector<T> hann(unsigned N, bool sym)
	{
		return general_hamming<T>(N, 0.5, sym);
	}

	template <class T>
	std::vector<T> bartlett(unsigned N, bool sym)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);

		auto n = arange<T>(0.0, static_cast<T>(M));
		auto w(n);

		for (auto& x : w)
		{
			if (x <= (M - 1) / 2.0)
			{
				x = static_cast<T>(2.0 * x / (M - 1));
			}
			else
			{
				x = static_cast<T>(2.0 - 2.0 * x / (M - 1));
			}
		}

		return Utilities::truncate(w, needs_trunc);
	}

	template<class T>
	std::vector<T> flattop(unsigned N, bool sym)
	{
		std::vector<T> a{
			static_cast<T>(0.21557895),
			static_cast<T>(0.41663158),
			static_cast<T>(0.277263158),
			static_cast<T>(0.083578947),
			static_cast<T>(0.006947368) };
		return general_cosine<T>(N, a, sym);
	}

	template <class T>
	std::vector<T> parzen(unsigned N, bool sym)
	{
		if (N <= 0) return std::vector<T>();

		auto [M, needs_trunc] = Utilities::extend(N, sym);
		auto n = arange<T>(-(static_cast<T>(M - 1)) / static_cast<T>(2.0),
			static_cast<T>((M - 1) / 2.0 + 0.5), 1.0);

		auto na = extract<T>([M](auto n) { return n < -static_cast<T>(M - 1) / 4.0; }, n);
		auto nb = extract<T>([M](auto n) {return std::abs(n) <= (M - 1) / 4.0; }, n);
		std::vector<T> wa;
		std::transform(na.begin(), na.end(), std::back_inserter(wa),
			[M](auto n) {return static_cast<T>(2*std::pow((1-std::abs(n) / static_cast<T>(M / 2.0)), 3)); });
		std::vector<T> wb;
		std::transform(nb.begin(), nb.end(), std::back_inserter(wb),
			[M](auto n)
			{
				return static_cast<T>(1 - 6 * std::pow(std::abs(n) / (M / 2.0), 2.0) + 6 * std::pow(std::abs(n) / (M / 2.0), 3.0));
			});

		auto w(wa);
		w.insert(w.end(), wb.begin(), wb.end());
		w.insert(w.end(), wa.rbegin(), wa.rend());
		
		return Utilities::truncate(w, needs_trunc);
		
	}


	// Explicit template instantiations
	template std::vector<float> boxcar(unsigned N, bool sym);
	template std::vector<double> boxcar(unsigned N, bool sym);
	template std::vector<long double> boxcar(unsigned N, bool sym);
	template std::vector<float> triang(unsigned N, bool sym);
	template std::vector<double> triang(unsigned N, bool sym);
	template std::vector<long double> triang(unsigned N, bool sym);
	template std::vector<float> general_cosine(unsigned N, const std::vector<float>& a, bool sym);
	template std::vector<double> general_cosine(unsigned N, const std::vector<double>& a, bool sym);
	template std::vector<long double> general_cosine(unsigned N, const std::vector<long double>& a, bool sym);
	template std::vector<float> general_hamming(unsigned N, double alpha, bool sym);
	template std::vector<double> general_hamming(unsigned N, double alpha, bool sym);
	template std::vector<long double> general_hamming(unsigned N, double alpha, bool sym);
	template std::vector<float> blackman(unsigned N, bool sym);
	template std::vector<double> blackman(unsigned N, bool sym);
	template std::vector<long double> blackman(unsigned N, bool sym);
	template std::vector<float> hamming(unsigned N, bool sym);
	template std::vector<double> hamming(unsigned N, bool sym);
	template std::vector<long double> hamming(unsigned N, bool sym);
	template std::vector<float> hann(unsigned N, bool sym);
	template std::vector<double> hann(unsigned N, bool sym);
	template std::vector<long double> hann(unsigned N, bool sym);
	template std::vector<float> bartlett(unsigned N, bool sym);
	template std::vector<double> bartlett(unsigned N, bool sym);
	template std::vector<long double> bartlett(unsigned N, bool sym);
	template std::vector<float> flattop(unsigned N, bool sym);
	template std::vector<double> flattop(unsigned N, bool sym);
	template std::vector<long double> flattop(unsigned N, bool sym);
	template std::vector<float> parzen(unsigned N, bool sym);
	template std::vector<double> parzen(unsigned N, bool sym);
	template std::vector<long double> parzen(unsigned N, bool sym);
}

