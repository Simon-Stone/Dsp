#include "Signal.h"

template<class T>
dsp::Signal<T>::Signal(const Signal& other) : samplingRate_Hz_(other.samplingRate_Hz_), samples_(other.samples_)
{
}

template <class T>
dsp::Signal<T>::Signal(Signal&& other) noexcept : samplingRate_Hz_(other.samplingRate_Hz_), samples_(std::move(other.samples_))
{
}

template <class T>
dsp::Signal<T>::Signal(unsigned samplingRate_Hz) : samplingRate_Hz_(samplingRate_Hz)
{
}

template <class T>
dsp::Signal<T>::Signal(const std::vector<T>& samples) : samples_(samples)
{
}

template <class T>
dsp::Signal<T>::Signal(unsigned samplingRate_Hz, const std::vector<T>& samples) : samplingRate_Hz_(samplingRate_Hz), samples_(samples)
{
}

template <class T>
const std::vector<T>& dsp::Signal<T>::getSamples() const
{
	return samples_;
}

template <class T>
std::vector<T>& dsp::Signal<T>::getSamples()
{
	return samples_;
}

template <class T>
unsigned dsp::Signal<T>::getSamplingRate_Hz() const
{
	return samplingRate_Hz_;
}

template<class T>
unsigned& dsp::Signal<T>::getSamplingRate_Hz()
{
	return samplingRate_Hz_;
}

template <class T>
void dsp::Signal<T>::setSamples(const std::vector<T>& samples)
{
	samples_ = samples;
}

template <class T>
void dsp::Signal<T>::setSamplingRate_Hz(unsigned newSamplingRate_Hz)
{
	samplingRate_Hz_ = newSamplingRate_Hz;
}

template <class T>
template <class U>
void dsp::Signal<T>::plus(U& x, const U& y)
{
	x += y;
}

template <class T>
template <class U>
void dsp::Signal<T>::plus(std::vector<U>& x, const U& y)
{
	std::transform(x.begin(), x.end(), x.begin(), [&y](auto& x) {return x + y; });
}

template <class T>
template <class U>
void dsp::Signal<T>::plus(std::vector<U>& x, const std::vector<U>& y)
{
	std::transform(y.begin(), y.end(), x.begin(), x.begin(), std::plus<U>());
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator+=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		plus(this->at(i), rhs[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator+=(const std::vector<T>& vec)
{
	if (this->size() != vec.size()) { throw std::logic_error("Signal and vector have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		plus(this->at(i), vec[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator+=(const_reference value)
{
	for (auto& x : *this)
	{
		plus(x, value);
	}	
	return *this;
}

template <class T>
template <class U>
void dsp::Signal<T>::minus(U& x, const U& y)
{
	x -= y;
}

template <class T>
template <class U>
void dsp::Signal<T>::minus(std::vector<U>& x, const U& y)
{
	std::transform(x.begin(), x.end(), x.begin(), [&y](auto& x) {return x - y; });
}

template <class T>
template <class U>
void dsp::Signal<T>::minus(std::vector<U>& x, const std::vector<U>& y)
{
	std::transform(y.begin(), y.end(), x.begin(), x.begin(), std::minus<>());
}


template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator-=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	for(size_type i = 0; i < this->size(); ++i)
	{
		minus(this->at(i), rhs[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator-=(const std::vector<T>& vec)
{
	if (this->size() != vec.size()) { throw std::logic_error("Signal and vector have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		minus(this->at(i), vec[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator-=(const_reference value)
{
	for (auto& x : *this)
	{
		minus(x, value);
	}
	return *this;
}

template <class T>
template <class U>
void dsp::Signal<T>::multiplies(U& x, const U& y)
{
	x *= y;
}

template <class T>
template <class U>
void dsp::Signal<T>::multiplies(std::vector<U>& x, const U& y)
{
	std::transform(x.begin(), x.end(), x.begin(), [&y](auto& x) {return x * y; });
}

template <class T>
template <class U>
void dsp::Signal<T>::multiplies(std::vector<U>& x, const std::vector<U>& y)
{
	std::transform(y.begin(), y.end(), x.begin(), x.begin(), std::multiplies<>());
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator*=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		multiplies(this->at(i), rhs[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator*=(const_reference value)
{
	for (auto& x : *this)
	{
		multiplies(x, value);
	}
	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator*=(const std::vector<T>& vec)
{
	if (this->size() != vec.size()) { throw std::logic_error("Signal and vector have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		multiplies(this->at(i), vec[i]);
	}

	return *this;
}

template <class T>
template <class U>
void dsp::Signal<T>::divides(U& x, const U& y)
{
	x /= y;
}

template <class T>
template <class U>
void dsp::Signal<T>::divides(std::vector<U>& x, const U& y)
{
	std::transform(x.begin(), x.end(), x.begin(), [&y](auto& x) {return x / y; });
}

template <class T>
template <class U>
void dsp::Signal<T>::divides(std::vector<U>& x, const std::vector<U>& y)
{
	std::transform(y.begin(), y.end(), x.begin(), x.begin(), std::divides<>());
}


template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator/=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		divides(this->at(i), rhs[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator/=(const std::vector<T>& vec)
{
	if (this->size() != vec.size()) { throw std::logic_error("Signal and vector have different lengths!"); }

	for (size_type i = 0; i < this->size(); ++i)
	{
		divides(this->at(i), vec[i]);
	}

	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator/=(const_reference value)
{
	for (auto& x : *this)
	{
		divides(x, value);
	}
	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::assign(size_type count, const T& value)
{
	samples_.assign(count, value); return *this;
}

template <class T>
typename dsp::Signal<T>::reference dsp::Signal<T>::at(size_type pos)
{
	return samples_.at(pos);
}

template <class T>
const T& dsp::Signal<T>::at(size_type pos) const
{
	return samples_.at(pos);
}

template <class T>
typename dsp::Signal<T>::value_type dsp::Signal<T>::getValue(size_type pos) const
{
	return samples_[pos];
}

template <class T>
typename dsp::Signal<T>::reference dsp::Signal<T>::front()
{
	return samples_.front();
}

template <class T>
const T& dsp::Signal<T>::front() const
{
	return samples_.front();
}

template <class T>
typename dsp::Signal<T>::reference dsp::Signal<T>::back()
{
	return samples_.back();
}

template <class T>
const T& dsp::Signal<T>::back() const
{
	return samples_.back();
}

template<class T>
T* dsp::Signal<T>::data()
{
	return samples_.data();
}

template <class T>
const T* dsp::Signal<T>::data() const
{
	return samples_.data();
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::begin()
{
	return samples_.begin();
}

template <class T>
typename dsp::Signal<T>::const_iterator dsp::Signal<T>::begin() const
{
	return samples_.begin();
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::end()
{
	return samples_.end();
}

template <class T>
typename dsp::Signal<T>::const_iterator dsp::Signal<T>::end() const
{
	return samples_.end();
}

template <class T>
typename dsp::Signal<T>::reverse_iterator dsp::Signal<T>::rbegin()
{
	return samples_.rbegin();
}

template <class T>
typename dsp::Signal<T>::const_reverse_iterator dsp::Signal<T>::rbegin() const
{
	return samples_.rbegin();
}

template <class T>
bool dsp::Signal<T>::empty() const
{
	return samples_.empty();
}

template <class T>
typename dsp::Signal<T>::size_type dsp::Signal<T>::size() const
{
	return samples_.size();
}

template <class T>
typename dsp::Signal<T>::size_type dsp::Signal<T>::max_size() const
{
	return samples_.max_size();
}

template <class T>
void dsp::Signal<T>::reserve(size_type new_cap)
{
	samples_.reserve(new_cap);
}

template <class T>
typename dsp::Signal<T>::size_type dsp::Signal<T>::capacity() const
{
	return samples_.capacity();
}

template <class T>
void dsp::Signal<T>::shrink_to_fit()
{
	samples_.shrink_to_fit();
}

template <class T>
void dsp::Signal<T>::clear()
{
	samples_.clear();
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::insert(const_iterator pos, const value_type& value)
{
	return samples_.insert(pos, value);
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::insert(const_iterator pos, value_type&& value)
{
	return samples_.insert(pos, value);
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::insert(const_iterator pos, size_type count, const value_type& value)
{
	return samples_.insert(pos, count, value);
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::insert(const_iterator pos, std::initializer_list<T> ilist)
{
	return samples_.insert(pos, ilist);
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::erase(const_iterator pos)
{
	return samples_.erase(pos);
}

template <class T>
typename dsp::Signal<T>::iterator dsp::Signal<T>::erase(const_iterator first, const_iterator last)
{
	return samples_.erase(first, last);
}

template <class T>
void dsp::Signal<T>::push_back(const T& value)
{
	samples_.push_back(value);
}

template <class T>
void dsp::Signal<T>::push_back(T&& value)
{
	samples_.push_back(value);
}

template <class T>
void dsp::Signal<T>::pop_back()
{
	samples_.pop_back();
}

template <class T>
void dsp::Signal<T>::resize(size_type count)
{
	samples_.resize(count);
}

template <class T>
void dsp::Signal<T>::resize(size_type count, const value_type value)
{
	samples_.resize(count, value);
}

template <class T>
void dsp::Signal<T>::swap(Signal& other) noexcept
{
	samples_.swap(other.getSamples());
}

template <class T>
typename dsp::Signal<T>::reference dsp::Signal<T>::operator[](size_type pos)
{
	return samples_[pos];
}

template <class T>
typename dsp::Signal<T>::const_reference dsp::Signal<T>::operator[](size_type pos) const
{
	return samples_[pos];
}

template <class T>
auto dsp::Signal<T>::real(const T& c)
{
	return std::real<T>(c);
}

template <class T>
template <class U>
auto dsp::Signal<T>::imag(const std::vector<U>& c)
{
	using V = decltype(imag(std::declval<U>()));
	std::vector<V> imagPart;

	std::transform(c.begin(), c.end(), std::back_inserter(imagPart.begin()), std::imag<V>);
	return imagPart;
}

template <class T>
auto dsp::Signal<T>::imag(const T& c)
{
	return std::imag<T>(c);
}

template <class T>
template <class U>
auto dsp::Signal<T>::real(const std::vector<U>& c)
{
	using V = decltype(real(std::declval<U>()));
	std::vector<V> realPart;

	std::transform(c.begin(), c.end(), std::back_inserter(realPart.begin()), std::real<V>);
	return realPart;
}


template <class T>
auto dsp::Signal<T>::real() const
{
	using U = decltype(real(std::declval<T>()));
	Signal<U> realPart(this->getSamplingRate_Hz());
	realPart.reserve(this->size());

	realPart.setSamplingRate_Hz(this->getSamplingRate_Hz());
	std::transform(samples_.begin(), samples_.end(), std::back_inserter(realPart), [this](auto c) {return Signal<T>::real(c); });
	return realPart;
}

template <class T>
auto dsp::Signal<T>::imag() const
{
	using U = decltype(real(std::declval<T>()));
	Signal<U> imagPart(this->getSamplingRate_Hz());
	imagPart.reserve(this->size());

	imagPart.setSamplingRate_Hz(this->getSamplingRate_Hz());
	std::transform(samples_.begin(), samples_.end(), std::back_inserter(imagPart), [this](auto c) {return Signal<T>::imag(c); });
	return imagPart;
}


template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator=(const Signal& other)
{
	if (this != &other)  // Self-assignment protection
	{
		samplingRate_Hz_ = other.getSamplingRate_Hz();
		samples_ = other.getSamples();
	}
	return *this;
}

template <class T>
dsp::Signal<T>& dsp::Signal<T>::operator=(Signal&& other) noexcept
{
	this->samplingRate_Hz_ = other.samplingRate_Hz_;
	this->samples_ = std::move(other.samples_);

	return *this;
}


// Integer type signals
template class dsp::Signal<short int>;
template class dsp::Signal<unsigned short int>;
template class dsp::Signal<int>;
template class dsp::Signal<unsigned int>;
template class dsp::Signal<long int>;
template class dsp::Signal<unsigned long int>;
template class dsp::Signal<long long int>;
template class dsp::Signal<unsigned long long int>;

// Floating point type signals
template class dsp::Signal<float>;
template class dsp::Signal<std::complex<float>>;
template class dsp::Signal<double>;
template class dsp::Signal<std::complex<double>>;
template class dsp::Signal<long double>;
template class dsp::Signal<std::complex<long double>>;

// Matrix-like signals
template class dsp::Signal<std::vector<float>>;
template class dsp::Signal<std::vector<double>>;
template class dsp::Signal<std::vector<long double>>;

