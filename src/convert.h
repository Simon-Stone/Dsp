#pragma once
#include <string>

/// @brief Conversions between various scales
namespace dsp::convert
{
	/// @brief Converts from Hertz to Bark scale
	/// @tparam T Data type of the values
	/// @param hz Frequency value in Hertz
	/// @return Frequency value in Bark
	template<class T>
	T hz2bark(T hz);

	/// @brief Converts from Bark to Hertz scale
	/// @tparam T Data type of the values
	/// @param bark Frequency value in Bark
	/// @return Frequency value in Hertz
	template<class T>
	T bark2hz(T bark);
	
	/// @brief Method to calculate the Mel frequencies
	enum class mel_method {slaney, stanley_smith, zwicker};
	
	/// @brief Converts from Hertz to Mel scale
	/// @tparam T Data type of the values
	/// @param hz Frequency value in Hertz
	/// @param method The method to calculate the Mel frequencies. 'stanley_smith' is the method used by HTK and Matlab (default). 'slaney' is the default method used by librosa. 'zwicker' is the linear conversion from Bark to Mel.
	/// @return Frequency value in Mel 
	template<class T>
	T hz2mel(T hz, mel_method method = mel_method::stanley_smith);

	/// @brief Converts from Mel to Hertz scale
	/// @tparam T Data type of the values
	/// @param method The method to calculate the Mel frequencies. 'stanley_smith' is the method used by HTK and Matlab (default). 'slaney' is the default method used by librosa. 'zwicker' is the linear conversion from Bark to Mel.
	/// @return Frequency value in Hertz
	template<class T>
	T mel2hz(T mel, mel_method method = mel_method::stanley_smith);

	/// @brief Converts from Hertz to MIDI note number
	/// @tparam T Data type of the values
	/// @param hz Frequency in Hertz
	/// @return MIDI note number
	template<class T>
	T hz2midi(T hz);

	/// @brief Converts from MIDI note number to Hertz
	/// @tparam T Data type of the values
	/// @param midi MIDI note number
	/// @return Frequency in Hertz
	template<class T>
	T midi2hz(T midi);

	/// @brief Converts a spelled note to the corresponding MIDI number
	/// @param note String spelling out the note. The following notations are allowed:
	/// Sharp: '#',
	/// Flat: 'b', '!',
	/// Natural: '',
	/// @param round_midi If true, round the exact number to nearest midi number.
	/// @return MIDI number
	double note2midi(std::string note, bool round_midi = true);

	/// @brief Converts a MIDI number to a note string
	/// @param midi MIDI number
	/// @param octave If true (default), include the octave number
	/// @param cents If true, cent markers will be appended for fractional notes, e.g., 'A4+03' Default: false.
	/// @param key key signature to use when resolving enharmonic equivalences
	/// @return Note string in the format 'C0', 'D#1', ...
	std::string midi2note(double midi, bool octave = true, bool cents = false, std::string key = "C:maj");
	
	/// @brief Frequency of the musical note C0 (octave 0) in Hertz
	static const double kFrequency_C0 = 16.35159783;

	/// @brief Converts from Hertz to semitones
	/// @tparam T Data type of the values
	/// @param hz Frequency in Hertz
	/// @param reference Reference frequency in Hertz. Default: 16.35159783 (musical note C0).
	/// @return Frequency value in semitones relative to reference
	template<class T>
	T hz2st(T hz, T reference = kFrequency_C0);

	/// @brief Converts from semitones to Hertz
	/// @tparam T Data type of the values
	/// @param st Frequency in semitones
	/// @param reference Reference frequency in Hertz. Default: 16.35159783 (musical note C0).
	/// @return Frequency value in Hertz
	template<class T>
	T st2hz(T st, T reference = kFrequency_C0);

}
