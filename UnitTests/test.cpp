#include "pch.h"

#include <numeric>
#include <iostream>
#include <chrono>

#include "../src/fft.h"
#include "../src/Signal.h"
#include "../src/window.h"
#include "../src/dsp.h"
#include "../src/filter.h"
#include "../src/signals.h"
#include "../src/convert.h"

struct DspTest : ::testing::Test
{

};

TEST_F(DspTest, ConvolutionTest)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000).getSamples();

	for (const auto& k : x)
	{
		std::cout << k << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	auto y = dsp::signals::sin<double>(100, 0.02, 8000, 1, dsp::pi / 2).getSamples();

	for (const auto& x : y)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	auto conv = dsp::convolve(x, y, dsp::convolution_mode::full);

	for (const auto& x : conv)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	conv = dsp::convolve(x, y, dsp::convolution_mode::valid);

	for (const auto& x : conv)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	conv = dsp::convolve(x, y, dsp::convolution_mode::same);

	for (const auto& x : conv)
	{
		std::cout << x << std::endl;
	}
}

TEST_F(DspTest, CorrelationTest)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000).getSamples();

	for (const auto& k : x)
	{
		std::cout << k << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	auto y = dsp::signals::cos<double>(100, 0.02, 8000).getSamples();

	for (const auto& x : y)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	auto corr = dsp::correlate(x, y, dsp::correlation_mode::full);

	for (const auto& x : corr)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	corr = dsp::correlate(x, y, dsp::correlation_mode::valid);

	for (const auto& x : corr)
	{
		std::cout << x << std::endl;
	}

	std::cout << "*******************************************" << std::endl;

	corr = dsp::correlate(x, y, dsp::correlation_mode::same);

	for (const auto& x : corr)
	{
		std::cout << x << std::endl;
	}
}

TEST_F(DspTest, ModuloTest)
{
	std::cout << dsp::mod(3.0, 5.0) << std::endl;

	std::cout << dsp::mod(1, 5) << std::endl;
}

TEST_F(DspTest, SignalEnergyTest)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000);

	for (const auto& k : x)
	{
		std::cout << k << std::endl;
	}

	auto e = dsp::calculateEnergy<double>(x.begin(), x.end());

	std::cout << e << std::endl;
}

TEST_F(DspTest, SignalPowerTest)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000);

	for (const auto& k : x)
	{
		std::cout << k << std::endl;
	}

	auto e = dsp::calculateMeanPower<double>(x.begin(), x.end());

	std::cout << e << std::endl;

}

TEST_F(DspTest, LpcTest)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000);

	auto a = dsp::filter::lpc(x.getSamples(), 12);

	for (const auto& ai : a)
	{
		std::cout << ai << std::endl;
	}
}

TEST_F(DspTest, FilterTest)
{
	std::vector<double> x{ 0.814723686393179, 0.905791937075619, 0.126986816293506,
		0.913375856139019, 0.632359246225410, 0.0975404049994095, 0.278498218867048,
		0.546881519204984, 0.957506835434298, 0.964888535199277, 0.157613081677548,
		0.970592781760616, 0.957166948242946, 0.485375648722841, 0.800280468888800 };

	auto y = dsp::filter::filter({ 1.0 }, { 1.0, -0.2 }, x);

	for (const auto& yi : y)
	{
		std::cout << yi << std::endl;
	}
}

TEST_F(DspTest, UnitConversions)
{
	auto f = 880.0;
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2bark(f);
	std::cout << f << " Bark" << std::endl;

	f = dsp::convert::bark2hz(f);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2mel(f);
	std::cout << f << " Mel (matlab)" << std::endl;

	f = dsp::convert::mel2hz(f);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2mel(f, dsp::convert::mel_method::slaney);
	std::cout << f << " Mel (librosa)" << std::endl;

	f = dsp::convert::mel2hz(f, dsp::convert::mel_method::slaney);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2mel(f, dsp::convert::mel_method::zwicker);
	std::cout << f << " Mel (zwicker)" << std::endl;

	f = dsp::convert::mel2hz(f, dsp::convert::mel_method::zwicker);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2midi(f);
	std::cout << f << " MIDI" << std::endl;

	f = dsp::convert::midi2hz(f);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2st(f);
	std::cout << f << " st (reference C0)" << std::endl;

	f = dsp::convert::st2hz(f);
	std::cout << f << " Hz" << std::endl;

	f = dsp::convert::hz2st(f, 440.0);
	std::cout << f << " st (reference A)" << std::endl;

	f = dsp::convert::st2hz(f);
	std::cout << f << " Hz" << std::endl;

	std::cout << dsp::convert::note2midi("C") << std::endl;
	std::cout << dsp::convert::note2midi("C#3") << std::endl;
	std::cout << dsp::convert::note2midi("f4") << std::endl;
	std::cout << dsp::convert::note2midi("Bb-1") << std::endl;
	std::cout << dsp::convert::note2midi("A!8") << std::endl;
}

TEST_F(DspTest, Window)
{
	auto win = dsp::window::get_window<double>(dsp::window::type::hann, 51);

	for (const auto& w : win)
	{
		std::cout << w << std::endl;
	}
}

TEST_F(DspTest, Spectrogram)
{
	auto x = dsp::signals::sin<double>(100, 0.02, 8000);

	auto X = dsp::fft::spectrogram(x.getSamples(), 256, 0.5, 8000);

	for (const auto& Xi : X)
	{
		for (const auto& Xik : Xi)
		{
			std::cout << Xik << std::endl;
		}
	}
}