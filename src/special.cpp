#include "special.h"

#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>

template<class T>
T dsp::special::i0(T x)
{
	return static_cast<T>(boost::math::cyl_bessel_i(0, x));
}


template<class T>
std::vector<T> dsp::special::i0(std::vector<T> x)
{
	std::vector<T> I;
	I.resize(x.size());

	std::transform(x.begin(), x.end(), I.begin(), [](auto x) {return i0(x); });
	return I;
}

// Explicit template instantiation
template float dsp::special::i0(float x);
template double dsp::special::i0(double x);
template long double dsp::special::i0(long double x);
template std::vector<float> dsp::special::i0<float>(std::vector<float> x);
template std::vector<double> dsp::special::i0<double>(std::vector<double> x);
template std::vector<long double> dsp::special::i0<long double>(std::vector<long double> x);