#pragma once
#include "Signal.h"
/// <summary>
/// An abstract base class for all kinds of transformation (e.g., Fast Fourier Transform)
/// </summary>

namespace Dsp
{
	template<class T>
	class Transform
	{
	public:
		Transform() = default;
		virtual ~Transform() = default;

		Signal<T> operator()(const Signal<T>& signal) { return transform(signal); }

	private:
		virtual Signal<T> transform(const Signal<T>& signal) = 0;
	};
}


