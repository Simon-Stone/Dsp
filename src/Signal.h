#pragma once
#include <algorithm>
#include <complex>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace Dsp
{
	/// <summary>
	/// This class represents a signal.
	///
	/// It is implemented as a facade design pattern wrapped around an STL vector.
	/// Fun fact: It is not derived from std::vector because that class does not have a virtual destructor, potentially producing memory leaks.
	/// </summary>
	///
	template <class T>
	class Signal
	{
	public:
		Signal() = default;
		Signal(const Signal& other);
		explicit Signal(unsigned samplingRate_Hz);
		explicit Signal(const std::vector<T>& samples);
		explicit Signal(unsigned samplingRate_Hz, const std::vector<T>& samples);
		virtual ~Signal() = default;

		// Getter
		const std::vector<T>& getSamples() const;
		std::vector<T>& getSamples();
		auto real() const;
		auto imag() const;
		unsigned getSamplingRate_Hz() const;

		// Setter
		void setSamples(const std::vector<T>& samples);
		void setSamplingRate_Hz(unsigned newSamplingRate_Hz);

		// Stream output
		friend std::ostream& operator<<(std::ostream& os, const Signal& obj)
		{
			for (const auto& x : obj.getSamples())
			{
				os << x << ",";
			}
			return os << std::endl;
		}

		/*
		 * Facade types
		*/
		using value_type = typename std::vector<T>::value_type;
		using allocator_type = typename std::vector<T>::allocator_type;
		using size_type = typename std::vector<T>::size_type;
		using difference_type = typename std::vector<T>::difference_type;
		using reference = typename std::vector<T>::reference;
		using const_reference = typename std::vector<T>::const_reference;
		using pointer = typename std::vector<T>::pointer;
		using const_pointer = typename std::vector<T>::const_pointer;
		using iterator = typename std::vector<T>::iterator;
		using const_iterator = typename std::vector<T>::const_iterator;
		using reverse_iterator = typename std::vector<T>::reverse_iterator;
		using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

		/*
		 * Arithmetic operations
		 */
		Signal<T>& operator+=(const Signal<T>& rhs);
		Signal<T>& operator+=(const_reference value);

		Signal<T>& operator-=(const Signal<T>& rhs);		
		Signal<T>& operator-=(const_reference value);

		Signal<T>& operator*=(const Signal<T>& rhs);
		Signal<T>& operator*=(const_reference value);

		Signal<T>& operator/=(const Signal<T>& rhs);
		Signal<T>& operator/=(const_reference value);
		
		// More operators are defined out of class for symmetry

	   /*
		* Facade methods
		*/
		
		Signal& assign(size_type count, const T& value);
		template< class InputIt >
		void assign(InputIt first, InputIt last);

		// Element access
		reference at(size_type pos);
		const T& at(size_type pos) const;

		reference front();
		const T& front() const;

		reference back();
		const T& back() const;

		T* data();
		const T* data() const;

		// Iterators
		iterator begin();
		const_iterator begin() const;

		iterator end();
		const_iterator end() const;

		reverse_iterator rbegin();
		const_reverse_iterator rbegin() const;

		// Capacity
		bool empty() const;
		size_type size() const;
		size_type max_size() const;
		void reserve(size_type new_cap);
		size_type capacity() const;
		void shrink_to_fit();

		// Modifiers
		void clear();

		iterator insert(const_iterator pos, const value_type& value);
		iterator insert(const_iterator pos, value_type&& value);
		iterator insert(const_iterator pos, size_type count, const value_type& value);
		template< class InputIt >
		iterator insert(const_iterator pos, InputIt first, InputIt last);
		iterator insert(const_iterator pos, std::initializer_list<T> ilist);

		template<class... Args>
		iterator emplace(const_iterator pos, Args&&... args);

		iterator erase(const_iterator pos);
		iterator erase(const_iterator first, const_iterator last);

		void push_back(const T& value);
		void push_back(T&& value);

		template< class... Args >
		void emplace_back(Args&&... args);

		void pop_back();

		void resize(size_type count);
		void resize(size_type count, const value_type value);

		void swap(Signal& other) noexcept;

		// Operators
		Signal& operator=(const Signal& other);
		reference& operator[](size_type pos);
		const_reference& operator[](size_type pos) const;
		// More operators are defined as non-member functions below		

	private:
		unsigned samplingRate_Hz_{ 0 };
		std::vector<T> samples_;  //!< A vector holding the actual samples		
	};
	
	template< class T>
	bool operator==(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		return (lhs.getSamplingRate_Hz() == rhs.getSamplingRate_Hz() && lhs.getSamples() == rhs.getSamples());
	}

	template< class T>
	bool operator!=(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		return !(lhs == rhs);
	}

	template< class T>
	bool operator<(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		if (lhs.getSamplingRate_Hz() != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
		return lhs.getSamples() < rhs.getSamples();
	}

	template <class T>
	bool operator<(const Signal<T>& lhs, const T& value)
	{
		// Look for the first element that is greater than or equal to value
		auto lb = std::lower_bound(lhs.begin(), lhs.end(), value);
		// If no such element was found, all elements are less than value
		return lb == lhs.end();
	}

	template< class T>
	bool operator<=(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		if (lhs.getSamplingRate_Hz() != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
		return lhs.getSamples() <= rhs.getSamples();
	}

	template <class T>
	bool operator<=(const Signal<T>& lhs, const T& value)
	{
		// Look for an element that is greater than value
		auto it = std::upper_bound(lhs.begin(), lhs.end(), value);
		// If no such element was found, all elements are less than or equal to value
		return it == lhs.end();
	}

	template< class T>
	bool operator>(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		if (lhs.getSamplingRate_Hz() != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
		return lhs.getSamples() > rhs.getSamples();
	}

	template <class T>
	bool operator>(const Signal<T>& lhs, const T& value)
	{
		// Look for an element that is smaller than or equal to value
		auto it = std::find_if(lhs.begin(), lhs.end(), [value](T x) {x <= value; });
		// If no such element was found, all elements are greater than value
		return it == lhs.end();
	}

	template< class T>
	bool operator>=(const Signal<T>& lhs, const Signal<T>& rhs)
	{
		if (lhs.getSamplingRate_Hz() != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
		return lhs.getSamples() >= rhs.getSamples();
	}

	template <class T>
	bool operator>=(const Signal<T>& lhs, const T& value)
	{
		// Look for an element that is smaller than value
		auto it = std::find_if(lhs.begin(), lhs.end(), [value](T x) {x < value; });
		// If no such element was found, all elements are greater than or equal to value
		return it == lhs.end();
	}

	template< class T>
	void swap(Signal<T>& lhs, Signal<T>& rhs) noexcept
	{
		if (lhs.getSamplingRate_Hz() != rhs.getSamplingRate_Hz()) { throw std::logic_error("Signals have different sampling rates!"); }
		std::vector<T>::swap(lhs.samples_, rhs.samples);
	}

	// Symmetric arithmetic operators
	template<class T>
	Signal<T> operator+(Signal<T>& lhs, const Signal<T>& rhs)
	{
		lhs += rhs;

		return lhs;
	}

	template<class T>
	Signal<T> operator+(Signal<T>& lhs, typename Signal<T>::const_reference value)
	{
		lhs += value;

		return lhs;
	}

	template<class T>
	Signal<T> operator-(Signal<T>& lhs, const Signal<T>& rhs)
	{
		lhs -= rhs;

		return lhs;
	}

	template<class T>
	Signal<T> operator-(Signal<T>& lhs, typename Signal<T>::const_reference value)
	{
		lhs -= value;

		return lhs;
	}

	template<class T>
	Signal<T> operator*(Signal<T>& lhs, const Signal<T>& rhs)
	{
		lhs *= rhs;

		return lhs;
	}

	template<class T>
	Signal<T> operator*(Signal<T>& lhs, typename Signal<T>::const_reference value)
	{
		lhs *= value;

		return lhs;
	}

	template<class T>
	Signal<T> operator/(Signal<T>& lhs, const Signal<T>& rhs)
	{
		lhs /= rhs;

		return lhs;
	}

	template<class T>
	Signal<T> operator/(Signal<T>& lhs, typename Signal<T>::const_reference value)
	{
		lhs /= value;

		return lhs;
	}

	template<class T>
	auto abs(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::abs() function applied to the signal's sample
		using U = decltype(std::abs(std::declval<T>()));
		Signal<U> absoluteSignal(signal.getSamplingRate_Hz());
		absoluteSignal.reserve(signal.size());
		
		std::transform(signal.begin(), signal.end(), std::back_inserter(absoluteSignal), [](auto x) {return std::abs(x); });

		return absoluteSignal;
	}

	template<class T>
	auto real(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::real() function applied to the signal's sample
		using U = decltype(std::real(std::declval<T>()));
		Signal<U> realPart(signal.getSamplingRate_Hz());
		realPart.reserve(signal.size());

		std::transform(signal.begin(), signal.end(), std::back_inserter(realPart), [](auto x) {return std::real(x); });

		return realPart;
	}

	template<class T>
	auto imag(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::imag() function applied to the signal's sample
		using U = decltype(std::imag(std::declval<T>()));
		Signal<U> imaginaryPart(signal.getSamplingRate_Hz());
		imaginaryPart.reserve(signal.size());

		std::transform(signal.begin(), signal.end(), std::back_inserter(imaginaryPart), [](auto x) {return std::imag(x); });

		return imaginaryPart;
	}

	template<class T>
	auto conj(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::conj() function applied to the signal's sample
		using U = decltype(std::conj(std::declval<T>()));
		Signal<U> conjugatedSignal(signal.getSamplingRate_Hz());
		conjugatedSignal.reserve(signal.size());

		std::transform(signal.begin(), signal.end(), std::back_inserter(conjugatedSignal), [](auto x) {return std::conj(x); });

		return conjugatedSignal;
	}

	template<class T>
	auto arg(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::arg() function applied to the signal's sample
		using U = decltype(std::arg(std::declval<T>()));
		Signal<U> phaseSignal(signal.getSamplingRate_Hz());
		phaseSignal.reserve(signal.size());

		std::transform(signal.begin(), signal.end(), std::back_inserter(phaseSignal), [](auto x) {return std::arg(x); });

		return phaseSignal;
	}

	template<class T>
	auto norm(Signal<T> signal)
	{
		// The returned signal should have the return type of the std::norm() function applied to the signal's sample
		using U = decltype(std::norm(std::declval<T>()));
		Signal<U> normSignal(signal.getSamplingRate_Hz());
		normSignal.reserve(signal.size());

		std::transform(signal.begin(), signal.end(), std::back_inserter(normSignal), [](auto x) {return std::norm(x); });

		return normSignal;
	}
}
