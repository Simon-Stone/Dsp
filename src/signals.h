#pragma once
#include "Signal.h"

/// @brief Convenience functions to generate various sampled test signals
namespace dsp::signals
{

	/// @brief Return a sampled sine signal
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param frequency_Hz The frequency of the sine in Hertz
	/// @param length_s The length of the signal in seconds
	/// @param samplingRate_Hz The sampling rate in Hertz
	/// @param amplitude Amplitude of the sine
	/// @param phase Starting phase angle of the sine in rad
	/// @return A signal containing the specified sine
	template <class T>
	Signal<T> sin(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude = 1.0, double phase = 0.0);

	/// @brief Return a sampled cosine signal
	/// @tparam T Type of the samples. Should be float, double, or long double. Other types will cause undefined behavior.
	/// @param frequency_Hz The frequency of the cosine in Hertz
	/// @param length_s The length of the signal in seconds
	/// @param samplingRate_Hz The sampling rate in Hertz
	/// @param amplitude Amplitude of the cosine
	/// @param phase Starting phase angle of the cosine in rad
	/// @return A signal containing the specified cosine
	template <class T>
	Signal<T> cos(unsigned frequency_Hz, double length_s, unsigned samplingRate_Hz, double amplitude = 1.0, double phase = 0.0);
}