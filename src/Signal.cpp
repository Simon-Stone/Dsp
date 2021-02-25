#include "Signal.h"

template<class T>
Dsp::Signal<T>::Signal(const Signal& other) : samplingRate_Hz_(other.getSamplingRate_Hz()), samples_(other.getSamples())
{
}

template <class T>
Dsp::Signal<T>::Signal(unsigned samplingRate_Hz) : samplingRate_Hz_(samplingRate_Hz)
{
}

template <class T>
Dsp::Signal<T>::Signal(const std::vector<T>& samples) : samples_(samples)
{
}

template <class T>
Dsp::Signal<T>::Signal(unsigned samplingRate_Hz, const std::vector<T>& samples) : samplingRate_Hz_(samplingRate_Hz), samples_(samples)
{	
}

template <class T>
const std::vector<T>& Dsp::Signal<T>::getSamples() const
{
	return samples_;
}

template <class T>
std::vector<T>& Dsp::Signal<T>::getSamples()
{
	return samples_;
}

template <class T>
unsigned Dsp::Signal<T>::getSamplingRate_Hz() const
{
	return samplingRate_Hz_;
}

template <class T>
void Dsp::Signal<T>::setSamples(const std::vector<T>& samples)
{
	samples_ = samples;
}

template <class T>
void Dsp::Signal<T>::setSamplingRate_Hz(unsigned newSamplingRate_Hz)
{
	samplingRate_Hz_ = newSamplingRate_Hz;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator+=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	std::transform(rhs.begin(), rhs.end(), this->begin(), this->begin(), std::plus<T>());

	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator+=(const_reference value)
{
	std::transform(this->begin(), this->end(), this->begin(), [&value](auto& y) {return y + value; });
	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator-=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	std::transform(rhs.begin(), rhs.end(), this->begin(), this->begin(), std::minus<T>());

	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator-=(const_reference value)
{
	std::transform(this->begin(), this->end(), this->begin(), [&value](auto& y) {return y - value; });
	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator*=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	std::transform(rhs.begin(), rhs.end(), this->begin(), this->begin(), std::multiplies<T>());

	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator*=(const_reference value)
{
	std::transform(this->begin(), this->end(), this->begin(), [&value](auto& y) {return y * value; });
	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator/=(const Signal<T>& rhs)
{
	if (this->samplingRate_Hz_ != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
	if (this->size() != rhs.size()) { throw std::logic_error("Signals have different lengths!"); }

	std::transform(rhs.begin(), rhs.end(), this->begin(), this->begin(), std::divides<T>());

	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator/=(const_reference value)
{
	std::transform(this->begin(), this->end(), this->begin(), [&value](auto& y) {return y / value; });
	return *this;
}

template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::assign(size_type count, const T& value)
{
	samples_.assign(count, value); return *this;
}

template <class T>
template <class InputIt>
void Dsp::Signal<T>::assign(InputIt first, InputIt last)
{
	samples_.assign(first, last);
}

template <class T>
typename Dsp::Signal<T>::reference Dsp::Signal<T>::at(size_type pos)
{
	return samples_.at(pos);
}

template <class T>
const T& Dsp::Signal<T>::at(size_type pos) const
{
	return samples_.at(pos);
}

template <class T>
typename Dsp::Signal<T>::reference Dsp::Signal<T>::front()
{
	return samples_.front();
}

template <class T>
const T& Dsp::Signal<T>::front() const
{
	return samples_.front();
}

template <class T>
typename Dsp::Signal<T>::reference Dsp::Signal<T>::back()
{
	return samples_.back();
}

template <class T>
const T& Dsp::Signal<T>::back() const
{
	return samples_.back();
}

template<class T>
T* Dsp::Signal<T>::data()
{
	return samples_.data();
}

template <class T>
const T* Dsp::Signal<T>::data() const
{
	return samples_.data();
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::begin()
{
	return samples_.begin();
}

template <class T>
typename Dsp::Signal<T>::const_iterator Dsp::Signal<T>::begin() const
{
	return samples_.begin();
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::end()
{
	return samples_.end();
}

template <class T>
typename Dsp::Signal<T>::const_iterator Dsp::Signal<T>::end() const
{
	return samples_.end();
}

template <class T>
typename Dsp::Signal<T>::reverse_iterator Dsp::Signal<T>::rbegin()
{
	return samples_.rbegin();
}

template <class T>
typename Dsp::Signal<T>::const_reverse_iterator Dsp::Signal<T>::rbegin() const
{
	return samples_.rbegin();
}

template <class T>
bool Dsp::Signal<T>::empty() const
{
	return samples_.empty();
}

template <class T>
typename Dsp::Signal<T>::size_type Dsp::Signal<T>::size() const
{
	return samples_.size();
}

template <class T>
typename Dsp::Signal<T>::size_type Dsp::Signal<T>::max_size() const
{
	return samples_.max_size();
}

template <class T>
void Dsp::Signal<T>::reserve(size_type new_cap)
{
	samples_.reserve(new_cap);
}

template <class T>
typename Dsp::Signal<T>::size_type Dsp::Signal<T>::capacity() const
{
	return samples_.capacity();
}

template <class T>
void Dsp::Signal<T>::shrink_to_fit()
{
	samples_.shrink_to_fit();
}

template <class T>
void Dsp::Signal<T>::clear()
{
	samples_.clear();
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::insert(const_iterator pos, const value_type& value)
{
	return samples_.insert(pos, value);
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::insert(const_iterator pos, value_type&& value)
{
	return samples_.insert(pos, value);
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::insert(const_iterator pos, size_type count, const value_type& value)
{
	return samples_.insert(pos, count, value);
}

template <class T>
template <class InputIt>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::insert(const_iterator pos, InputIt first, InputIt last)
{
	return samples_.insert(pos, first, last);
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::insert(const_iterator pos, std::initializer_list<T> ilist)
{
	return samples_.insert(pos, ilist);
}

template <class T>
template <class ... Args>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::emplace(const_iterator pos, Args&&... args)
{
	return samples_.emplace(pos, args);
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::erase(const_iterator pos)
{
	return samples_.erase(pos);
}

template <class T>
typename Dsp::Signal<T>::iterator Dsp::Signal<T>::erase(const_iterator first, const_iterator last)
{
	return samples_.erase(first, last);
}

template <class T>
void Dsp::Signal<T>::push_back(const T& value)
{
	samples_.push_back(value);
}

template <class T>
void Dsp::Signal<T>::push_back(T&& value)
{
	samples_.push_back(value);
}

template <class T>
template <class ... Args>
void Dsp::Signal<T>::emplace_back(Args&&... args)
{
	samples_.emplace_back(args);
}

template <class T>
void Dsp::Signal<T>::pop_back()
{
	samples_.pop_back();
}

template <class T>
void Dsp::Signal<T>::resize(size_type count)
{
	samples_.resize(count);
}

template <class T>
void Dsp::Signal<T>::resize(size_type count, const value_type value)
{
	samples_.resize(count, value);
}

template <class T>
void Dsp::Signal<T>::swap(Signal& other) noexcept
{
	samples_.swap(other.getSamples());
}

template <class T>
typename Dsp::Signal<T>::reference& Dsp::Signal<T>::operator[](size_type pos)
{
	return samples_[pos];
}

template <class T>
typename Dsp::Signal<T>::const_reference& Dsp::Signal<T>::operator[](size_type pos) const
{
	return samples_[pos];
}

template <class T>
auto Dsp::Signal<T>::real() const
{
	// The returned signal should have the return type of the std::real() function applied to the called signal's sample
	using U = decltype(std::real(std::declval<T>()));
	Signal<U> realPart(this->getSamplingRate_Hz());
	realPart.reserve(this->size());

	realPart.setSamplingRate_Hz(this->getSamplingRate_Hz());
	std::transform(samples_.begin(), samples_.end(), std::back_inserter(realPart), [](auto c) {return std::real(c); });

	return realPart;
}

template <class T>
auto Dsp::Signal<T>::imag() const
{
	// The returned signal should have the return type of the std::imag() function applied to the called signal's sample
	using U = decltype(std::imag(std::declval<T>()));
	Dsp::Signal<U> imaginaryPart(this->getSamplingRate_Hz());
	imaginaryPart.reserve(this->size());

	imaginaryPart.setSamplingRate_Hz(this->getSamplingRate_Hz());
	std::transform(samples_.begin(), samples_.end(), std::back_inserter(imaginaryPart), [](auto c) {return std::imag(c); });

	return imaginaryPart;
}


template <class T>
Dsp::Signal<T>& Dsp::Signal<T>::operator=(const Signal& other)
{
	if (this != &other)  // Self-assignment protection
	{
		samplingRate_Hz_ = other.getSamplingRate_Hz();
		samples_ = other.getSamples();
	}
	return *this;
}


// Integer type signals
template class Dsp::Signal<short int>;
template class Dsp::Signal<unsigned short int>;
template class Dsp::Signal<int>;
template class Dsp::Signal<unsigned int>;
template class Dsp::Signal<long int>;
template class Dsp::Signal<unsigned long int>;
template class Dsp::Signal<long long int>;
template class Dsp::Signal<unsigned long long int>;

// Floating point type signals
template class Dsp::Signal<float>;
template class Dsp::Signal<std::complex<float>>;
template class Dsp::Signal<double>;
template class Dsp::Signal<std::complex<double>>;
template class Dsp::Signal<long double>;
template class Dsp::Signal<std::complex<long double>>;

