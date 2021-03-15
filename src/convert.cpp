#include "convert.h"

#include <cmath>
#include <iostream>
#include <map>
#include <regex>
#include <stdexcept>

template <class T>
T dsp::convert::hz2bark(T hz)
{
	auto bark = 26.81 * hz / (1960.0 + hz) - 0.53;
	if (bark < 2.0)
	{
		bark = bark + 0.15 * (2 - bark);
	}
	if (bark > 20.1)
	{
		bark = bark + 0.22 * (bark - 20.1);
	}

	return static_cast<T>(bark);
}

template <class T>
T dsp::convert::bark2hz(T bark)
{
	if (bark < 2)
	{
		bark = static_cast<T>((bark - 0.3) / 0.85);
	}
	if (bark > 20.1)
	{
		bark = static_cast<T>((bark + 4.422) / 1.22);
	}

	auto hz = 1960.0 * (bark + 0.53) / (26.28 - bark);
	return static_cast<T>(hz);
}

template <class T>
T dsp::convert::hz2mel(T hz, mel_method method)
{
	T mel;
	switch (method)
	{
	case mel_method::slaney:
	{
		/* Slaney composes the hz2mel curve from a linear part and a logarithmic part */
		const T f_min{ 0.0 };
		const T f_slope{ static_cast<T>(200.0 / 3.0) };
		mel = (hz - f_min) / f_slope;

		T min_log_hz{ 1000.0 };  // Start of the log area in Hz
		T min_log_mel{ (min_log_hz - f_min) / f_slope };  // Start of the log area in Mel
		T logstep{ static_cast<T>(std::log(6.4) / 27.0) };

		if (hz >= min_log_hz)
		{
			mel = min_log_mel + std::log(hz / min_log_hz) / logstep;
		}
	}
	break;
	case mel_method::stanley_smith:
		mel = static_cast<T>(2595.0 * std::log10(1.0 + hz / 700.0));
		break;
	case mel_method::zwicker:
		mel = static_cast<T>(hz2bark(hz) * 100.0);
		break;
	default:
		throw std::runtime_error("Unknown method to calculate Mel frequencies!");
	}
	return static_cast<T>(mel);
}

template <class T>
T dsp::convert::mel2hz(T mel, mel_method method)
{
	T hz;
	switch (method)
	{
	case mel_method::slaney:
	{
		/* Slaney composes the mel2hz curve from a linear part and a logarithmic part */
		T f_min{ 0.0 };
		T f_slope{ static_cast<T>(200.0 / 3.0) };
		hz = f_min + f_slope * mel;

		T min_log_hz{ 1000.0 };  // Start of the log area in Hz
		T min_log_mel{ (min_log_hz - f_min) / f_slope };  // Start of the log area in Mel
		T logstep{ static_cast<T>(std::log(6.4) / 27.0) };

		if (mel >= min_log_mel)
		{
			hz = static_cast<T>(min_log_hz * std::exp(logstep * (mel - min_log_mel)));
		}
	}
	break;
	case mel_method::stanley_smith:
		hz = static_cast<T>(700.0 * (std::pow(10, mel / 2595.0) - 1.0));
		break;
	case mel_method::zwicker:
		hz = static_cast<T>(bark2hz(mel / 100.0));
		break;
	default:
		throw std::runtime_error("Unknown method to calculate Mel frequencies!");
	}

	return static_cast<T>(hz);
}

template <class T>
T dsp::convert::hz2midi(T hz)
{
	return static_cast<T>(12.0 * (std::log2(hz) - std::log2(440.0)) + 69.0);
}

template <class T>
T dsp::convert::midi2hz(T midi)
{
	return static_cast<T>(440.0 * std::pow(2.0, (midi - 69.0) / 12.0));
}

double dsp::convert::note2midi(std::string note, bool round_midi)
{
	std::map<std::string, int> pitch_map{ {"C", 0}, {"D", 2}, {"E", 4},
		{"F", 5}, {"G", 7}, {"A", 9}, {"B", 11} };
	std::map<std::string, int> acc_map{
		{"#", 1},
		{"", 0},
		{"b", -1},
		{"!", -1},
	};

	std::regex r(u8R"(^([A-Ga-g])([#♯𝄪b!♭𝄫♮]*)([+-]?[0-9]+)([+-][0-9]+)?$)");
	std::smatch m;
	std::regex_search(note, m, r);

	if (m.size() == 1) { throw std::runtime_error("Badly formatted note!"); }

	auto pitch = m[1].str();
	for (auto& c : pitch)
	{
		c = std::toupper(c);
	}
	auto offset = 0;
	for (auto& o : m[2].str())
	{
		auto acc = acc_map.find(std::string(1, o));
		if (acc == acc_map.end()) { throw std::runtime_error("Unknown accidental symbol!"); }
		offset += acc->second;
	}
	int octave = 0;
	if (m.size() > 3)
	{
		octave = std::stoi(m[3].str());
	}
	double cents = 0;
	if(m.size() > 3)
	{
		cents = std::stoi(m[3].str()) * 1e-2;
	}
	auto note_value = 12.0 * (octave + 1) + pitch_map[pitch] + offset + cents;

	if (round_midi)
	{
		note_value = std::round(note_value);
	}
	
	return note_value;
}

std::string dsp::convert::midi2note(double midi, bool octave, bool cents, std::string key)
{
	if (cents && !octave) { throw std::runtime_error("Cannot encode cents without octave information!"); }

	// TODO: Re-implement midi_to_note() from librosa: https://librosa.org/doc/0.8.0/_modules/librosa/core/convert.html#midi_to_note
	throw std::runtime_error("midi2note() not implemented yet!");

	return "";
}

template <class T>
T dsp::convert::hz2st(T hz, T reference)
{
	if (hz < 1.0)
	{
		hz = 1.0;
	}

	return static_cast<T>(12.0 * std::log(hz / reference) / std::log(2.0));
}

template <class T>
T dsp::convert::st2hz(T st, T reference)
{
	return static_cast<T>(reference * std::pow(2, st / 12.0));
}


// Explicit template instantiations
template float dsp::convert::hz2bark(float hz);
template double dsp::convert::hz2bark(double hz);
template long double dsp::convert::hz2bark(long double hz);

template float dsp::convert::bark2hz(float bark);
template double dsp::convert::bark2hz(double bark);
template long double dsp::convert::bark2hz(long double bark);

template float dsp::convert::hz2mel(float hz, mel_method method);
template double dsp::convert::hz2mel(double hz, mel_method method);
template long double dsp::convert::hz2mel(long double hz, mel_method method);

template float dsp::convert::mel2hz(float mel, mel_method method);
template double dsp::convert::mel2hz(double mel, mel_method method);
template long double dsp::convert::mel2hz(long double mel, mel_method method);

template float dsp::convert::hz2midi(float hz);
template double dsp::convert::hz2midi(double hz);
template long double dsp::convert::hz2midi(long double hz);

template float dsp::convert::midi2hz(float midi);
template double dsp::convert::midi2hz(double midi);
template long double dsp::convert::midi2hz(long double midi);

template float dsp::convert::hz2st(float hz, float reference);
template double dsp::convert::hz2st(double hz, double reference);
template long double dsp::convert::hz2st(long double hz, long double reference);

template float dsp::convert::st2hz(float st, float reference);
template double dsp::convert::st2hz(double st, double reference);
template long double dsp::convert::st2hz(long double st, long double reference);