#include "gtest/gtest.h"

#include <numeric>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>


#include "dsp.h"

struct DspTest : ::testing::Test
{

};

TEST_F(DspTest, SignalOperations)
{
	auto c = dsp::signals::cos<double>(100, 0.02, 8000);
	auto s = dsp::signals::sin<double>(100, 0.02, 8000);

	// 1 = cos^2(x) + sin^2(x) 
	auto sum = dsp::pow(c, 2) + dsp::pow(s, 2);
	auto ones = dsp::Signal<double>(sum.getSamplingRate_Hz(), dsp::signals::ones<double>(sum.size()));
	EXPECT_TRUE((sum - ones) < 0.00001);

	// cos(2x) = cos^2(x) - sin^2(x)
	auto diff = dsp::pow(c, 2) - dsp::pow(s, 2);
	auto c2 = dsp::signals::cos<double>(200, 0.02, 8000);

	EXPECT_TRUE((diff - c2) < 0.00001);
}

TEST_F(DspTest, Mean)
{
	auto c = dsp::signals::cos<double>(100, 0.02, 8000);

	EXPECT_NEAR(dsp::mean<double>(c.begin(), c.end()), 0.0, 1e-16);
}

TEST_F(DspTest, Median)
{
	std::vector<int> v{ 3, 2, 2, 1, 7 };
	std::vector<double> v2{ 3, 2, 3, 1 };
	auto s = dsp::Signal(v);
	auto s2 = dsp::Signal(v2);
	// Iterator interface
	EXPECT_EQ(dsp::median<int>(v.begin(), v.end()), 2);
	EXPECT_EQ(dsp::median<double>(v2.begin(), v2.end()), 2.5);
	EXPECT_EQ(dsp::median<int>(s.begin(), s.end()), 2);
	EXPECT_EQ(dsp::median<double>(s2.begin(), s2.end()), 2.5);

	// Reference interface
	EXPECT_EQ(dsp::median(v), 2);
	EXPECT_EQ(dsp::median(s), 2);
	EXPECT_EQ(dsp::median(v2), 2.5);
	EXPECT_EQ(dsp::median(s2), 2.5);
}

TEST_F(DspTest, Mode)
{
	std::vector<int> v{ 3, 2, 2, 1, 7 };

	EXPECT_EQ(dsp::mode<int>(v.begin(), v.end()), 2);
	EXPECT_EQ(dsp::mode(v), 2);
	EXPECT_EQ(dsp::mode(dsp::Signal<int>(v)), 2);
}

TEST_F(DspTest, Pad)
{
	std::vector<int> v{ 1,2,4,5 };
	auto v_pad = dsp::pad(v, { 2, 3 });

	for (const auto& x : v_pad)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;

	v_pad = dsp::pad(v, { 3, 4 }, { 2, 8 });

	for (const auto& x : v_pad)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;

	dsp::Signal<int> s(100, { 2, 3, 4, 5 });

	auto s_pad = dsp::pad(s, { 1, 2 }, { 1, 0 });

	for (const auto& x : s_pad)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;
}

TEST_F(DspTest, MedianFilter)
{
	std::vector<double> v{ 2.0, 3.0, 3.0, 4.0, 5.0, 3.0, 2.0 };
	auto v_med = dsp::filter::medianfilter(v, 5);

	for (const auto& x : v_med)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;

	EXPECT_EQ(v.size(), v_med.size());

	dsp::Signal<double> s(100, { 2.0, 3.0, 3.0, 4.0, 5.0, 3.0, 2.0 });
	auto s_med = dsp::filter::medianfilter(s, 5);

	for (const auto& x : s_med)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;

	EXPECT_EQ(s.size(), s_med.size());

}

TEST_F(DspTest, Unique)
{
	std::vector<int> v{ 3, 4, 2, 2, 1, 7 };

	auto u = dsp::unique(dsp::Signal(v));

	for (const auto& x : u)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;

	u = dsp::unique<int>(v.begin(), v.end());
	for (const auto& x : u)
	{
		std::cout << x << " ";
	}
	std::cout << std::endl;
}

TEST_F(DspTest, Zscore)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(5.0, 2.0);
	const int nsamples = 100000;
	dsp::Signal<double> x;
	for (int i = 0; i < nsamples; ++i)
	{
		x.push_back(distribution(generator));
	}

	auto standardized_x_vec = dsp::zscore<double>(x.begin(), x.end());
	EXPECT_NEAR(dsp::mean<double>(standardized_x_vec), 0, 0.000001);
	EXPECT_NEAR(dsp::std<double>(standardized_x_vec), 1, 0.000001);

	auto standardized_x_vec2 = dsp::zscore(x.getSamples());
	EXPECT_NEAR(dsp::mean<double>(standardized_x_vec2), 0, 0.000001);
	EXPECT_NEAR(dsp::std<double>(standardized_x_vec2), 1, 0.000001);

	auto standardized_x_sig = dsp::zscore(x);
	EXPECT_NEAR(dsp::mean<double>(standardized_x_vec2), 0, 0.000001);
	EXPECT_NEAR(dsp::std<double>(standardized_x_vec2), 1, 0.000001);
	EXPECT_EQ(standardized_x_sig.getSamplingRate_Hz(), x.getSamplingRate_Hz());
}

TEST_F(DspTest, VarAndStd)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(5.0, 2.0);
	const int nsamples = 100000;
	std::vector<double> x;
	for (int i = 0; i < nsamples; ++i)
	{
		x.push_back(distribution(generator));
	}

	auto mean = dsp::mean(x);
	EXPECT_NEAR(mean, 5.0, 0.1);
	auto var = dsp::var(x);
	EXPECT_NEAR(var, 4.0, 0.1);
}

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
	auto win = dsp::window::get_window<double>(dsp::window::type::gaussian, 128, true, { 32.0 });

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

TEST_F(DspTest, Benchmarking)
{
	std::ofstream outFile("../scripts/dsp_benchmark.csv", std::ios::out);

	outFile << "Task\tImplementation\tExecution time [us]" << std::endl;

	auto x = 7 * dsp::signals::sin<double>(100, 1, 48000);

	std::cout << "*********************************************************" << std::endl;
	std::cout << "********             z-score                  ***********" << std::endl;
	std::cout << "*********************************************************" << std::endl;
	for (int j = 0; j < 100; ++j)
	{
		// Naive implementation
		auto start = std::chrono::high_resolution_clock::now();
		double sum{ 0.0 };
		for (int i = 0; i < x.size(); ++i)
		{
			sum += x[i];
		}
		double x_mean = sum / static_cast<double>(x.size());

		sum = 0.0;
		for (int i = 0; i < x.size(); ++i)
		{
			sum += std::pow(x[i] - x_mean, 2);
		}
		double std = std::sqrt(sum) / (x.size() - 1);

		std::vector<double> standardized_x;
		for (int i = 0; i < x.size(); ++i)
		{
			standardized_x.push_back((x[i] - x_mean) / std);
		}

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		std::cout << "Naive implementation: " << "\t\t" << duration.count() << " µs" << std::endl;

		outFile << "zscore\tnaive\t" << duration.count() << std::endl;

		// DSP lib implementation
		auto start_dsp = std::chrono::high_resolution_clock::now();
		auto standardized_x_dsp = dsp::zscore(x);

		auto stop_dsp = std::chrono::high_resolution_clock::now();
		auto duration_dsp = std::chrono::duration_cast<std::chrono::microseconds>(stop_dsp - start_dsp);
		std::cout << "Modern DSP implementation: " << "\t" << duration_dsp.count() << " µs" << std::endl;

		outFile << "zscore\tdsp\t" << duration_dsp.count() << std::endl;
	}

	std::cout << "*********************************************************" << std::endl;
	std::cout << "********             energy                   ***********" << std::endl;
	std::cout << "*********************************************************" << std::endl;
	for (int j = 0; j < 100; ++j)
	{
		// Naive implementation
		auto start = std::chrono::high_resolution_clock::now();
		double sum{ 0.0 };
		for (int i = 0; i < x.size(); ++i)
		{
			sum += x[i] * x[i];
		}
		double x_energy = sum / static_cast<double>(x.size());
		auto stop = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		std::cout << "Naive implementation: " << "\t\t" << duration.count() << " µs" << std::endl;

		outFile << "energy\tnaive\t" << duration.count() << std::endl;
		// DSP lib implementation
		auto start_dsp = std::chrono::high_resolution_clock::now();
		auto x_energy_dsp = dsp::calculateMeanPower<double>(x.begin(), x.end());

		auto stop_dsp = std::chrono::high_resolution_clock::now();
		auto duration_dsp = std::chrono::duration_cast<std::chrono::microseconds>(stop_dsp - start_dsp);
		std::cout << "Modern DSP implementation: " << "\t" << duration_dsp.count() << " µs" << std::endl;

		outFile << "energy\tdsp\t" << duration_dsp.count() << std::endl;
	}

	std::cout << "*********************************************************" << std::endl;
	std::cout << "********   log-squared magnitude spectrum     ***********" << std::endl;
	std::cout << "*********************************************************" << std::endl;
	for (int j = 0; j < 100; ++j)
	{
		// Naive implementation
		auto start = std::chrono::high_resolution_clock::now();

		auto X = dsp::fft::rfft<double>({ x.begin(), x.end() }, x.size(), dsp::fft::NormalizationMode::backward, dsp::fft::backend::simple);

		std::vector<double> magSpec;
		for (int i = 0; i < X.size(); ++i)
		{
			magSpec.push_back(10 * std::log10(std::abs(X[i]) * std::abs(X[i])));
		}

		auto stop = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		std::cout << "Naive implementation: " << "\t\t" << duration.count() << " µs" << std::endl;
		outFile << "magspec\tnaive\t" << duration.count() << std::endl;

		// DSP lib implementation
		auto start_dsp = std::chrono::high_resolution_clock::now();
		auto x_magSpec = dsp::fft::logSquaredMagnitudeSpectrum(x.getSamples(), x.size(), 1);

		auto stop_dsp = std::chrono::high_resolution_clock::now();
		auto duration_dsp = std::chrono::duration_cast<std::chrono::microseconds>(stop_dsp - start_dsp);
		std::cout << "Modern DSP implementation: " << "\t" << duration_dsp.count() << " µs" << std::endl;

		outFile << "magspec\tdsp\t" << duration_dsp.count() << std::endl;
	}

}


TEST_F(DspTest, DctTest)
{
	// Define test vectors for each data type.
	const std::vector<int> VInt1{ 4, 3 };

	const unsigned int nBasisVectors = 20;
	double epsi = 1e-5;

	// Define test input - output pairs.
	const std::vector<long double> xLongDouble{0.814723686393179,
										   0.905791937075619,
										   0.126986816293506,
										   0.913375856139019,
										   0.632359246225410,
										   0.097540404999410,
										   0.278498218867048,
										   0.546881519204984,
										   0.957506835434298,
										   0.964888535199277,
										   0.157613081677548,
										   0.970592781760616,
										   0.957166948242946,
										   0.485375648722841,
										   0.800280468888800,
										   0.141886338627215,
										   0.421761282626275,
										   0.915735525189067,
										   0.792207329559554,
										   0.959492426392903};

	const std::vector<long double> yLongDouble{2.871259956478833,
											-0.199751167348291,
											 0.111463876686767,
											 0.111767508951061,
											 0.652160412955051,
											-0.186013450897578,
											 0.114716368511974,
											-0.231602462563853,
											-0.195967297690951,
											 0.625866739142916,
											 0.400441527850273,
											 0.092873866932537,
											-0.231641556908577,
											 0.030554201728932,
											 0.125354091788687,
											-0.831943488191366,
											-0.212855574804055,
											 0.025515454568773,
											-0.237438360840790,
											-0.113376431328529};


	const std::vector<double> xDouble(xLongDouble.begin(), xLongDouble.end());
	const std::vector<double> yDouble(yLongDouble.begin(), yLongDouble.end());

	const std::vector<float> xFloat(xLongDouble.begin(), xLongDouble.end());
	const std::vector<float> yFloat(yLongDouble.begin(), yLongDouble.end());

	// Define test target(s) for the cosineBasisVector calculation.
	std::vector<std::vector<double>> basisVectors =
	{
		{0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979, 0.223606797749979},
		{0.315252941349890, 0.307490367669328, 0.292156360634725, 0.269628494389924, 0.240461479767574, 0.205373505467449, 0.165228553867129, 0.121015126908468, 0.0738219058991409, 0.0248109445657177, -0.0248109445657176, -0.0738219058991408, -0.121015126908468, -0.165228553867129, -0.205373505467449, -0.240461479767574, -0.269628494389924, -0.292156360634725, -0.307490367669328, -0.315252941349890},
		{0.312334477467278, 0.281761002650515, 0.223606797749979, 0.143564401525505, 0.0494689214077114, -0.0494689214077113, -0.143564401525505, -0.223606797749979, -0.281761002650515, -0.312334477467278, -0.312334477467278, -0.281761002650515, -0.223606797749979, -0.143564401525505, -0.0494689214077114, 0.0494689214077113, 0.143564401525505, 0.223606797749979, 0.281761002650515, 0.312334477467278},
		{0.307490367669328, 0.240461479767574, 0.121015126908468, -0.0248109445657176, -0.165228553867129, -0.269628494389924, -0.315252941349890, -0.292156360634725, -0.205373505467449, -0.0738219058991411, 0.0738219058991409, 0.205373505467449, 0.292156360634725, 0.315252941349890, 0.269628494389924, 0.165228553867129, 0.0248109445657176, -0.121015126908468, -0.240461479767574, -0.307490367669328},
		{0.300750477503773, 0.185874017230092, 1.93633660727019e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -5.80900982181058e-17, 0.185874017230092, 0.300750477503773, 0.300750477503773, 0.185874017230092, 9.68168303635097e-17, -0.185874017230092, -0.300750477503773, -0.300750477503773, -0.185874017230092, -1.35543562508914e-16, 0.185874017230092, 0.300750477503773},
		{0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908468, -0.292156360634725, -0.292156360634725, -0.121015126908468, 0.121015126908468, 0.292156360634725, 0.292156360634725, 0.121015126908468, -0.121015126908467, -0.292156360634725},
		{0.281761002650515, 0.0494689214077114, -0.223606797749979, -0.312334477467278, -0.143564401525505, 0.143564401525504, 0.312334477467278, 0.223606797749979, -0.0494689214077115, -0.281761002650515, -0.281761002650515, -0.0494689214077112, 0.223606797749979, 0.312334477467278, 0.143564401525505, -0.143564401525504, -0.312334477467278, -0.223606797749979, 0.0494689214077109, 0.281761002650515},
		{0.269628494389924, -0.0248109445657176, -0.292156360634725, -0.240461479767574, 0.0738219058991409, 0.307490367669328, 0.205373505467449, -0.121015126908468, -0.315252941349890, -0.165228553867129, 0.165228553867129, 0.315252941349890, 0.121015126908468, -0.205373505467449, -0.307490367669328, -0.0738219058991409, 0.240461479767574, 0.292156360634725, 0.0248109445657175, -0.269628494389924},
		{0.255833636800846, -0.0977197537924274, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924273, -0.316227766016838, -0.0977197537924275, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924276, 0.255833636800846, 0.255833636800847, -0.0977197537924272, -0.316227766016838, -0.0977197537924277, 0.255833636800846},
		{0.240461479767574, -0.165228553867129, -0.292156360634725, 0.0738219058991409, 0.315252941349890, 0.0248109445657179, -0.307490367669328, -0.121015126908468, 0.269628494389924, 0.205373505467450, -0.205373505467449, -0.269628494389924, 0.121015126908468, 0.307490367669328, -0.0248109445657181, -0.315252941349890, -0.0738219058991410, 0.292156360634724, 0.165228553867129, -0.240461479767574},
		{0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749979, -0.223606797749979, -0.223606797749979, 0.223606797749979, 0.223606797749978, -0.223606797749979, -0.223606797749980, 0.223606797749979},
		{0.205373505467449, -0.269628494389924, -0.121015126908468, 0.307490367669328, 0.0248109445657176, -0.315252941349890, 0.0738219058991411, 0.292156360634725, -0.165228553867129, -0.240461479767575, 0.240461479767574, 0.165228553867129, -0.292156360634725, -0.0738219058991421, 0.315252941349890, -0.0248109445657169, -0.307490367669328, 0.121015126908467, 0.269628494389924, -0.205373505467448},
		{0.185874017230092, -0.300750477503773, -5.80900982181058e-17, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, 7.36003649626590e-16, -0.300750477503773, 0.185874017230092, 0.185874017230092, -0.300750477503773, -8.52183846062801e-16, 0.300750477503773, -0.185874017230092, -0.185874017230093, 0.300750477503773, -1.55102667445532e-16, -0.300750477503773, 0.185874017230092},
		{0.165228553867129, -0.315252941349890, 0.121015126908468, 0.205373505467450, -0.307490367669328, 0.0738219058991406, 0.240461479767574, -0.292156360634725, 0.0248109445657181, 0.269628494389924, -0.269628494389924, -0.0248109445657175, 0.292156360634725, -0.240461479767574, -0.0738219058991411, 0.307490367669328, -0.205373505467450, -0.121015126908469, 0.315252941349890, -0.165228553867129},
		{0.143564401525505, -0.312334477467278, 0.223606797749979, 0.0494689214077117, -0.281761002650515, 0.281761002650515, -0.0494689214077120, -0.223606797749980, 0.312334477467278, -0.143564401525504, -0.143564401525504, 0.312334477467278, -0.223606797749979, -0.0494689214077125, 0.281761002650516, -0.281761002650515, 0.0494689214077106, 0.223606797749980, -0.312334477467278, 0.143564401525503},
		{0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468, -0.121015126908468, 0.292156360634725, -0.292156360634725, 0.121015126908468, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908467, -0.121015126908469, 0.292156360634725, -0.292156360634724, 0.121015126908467, 0.121015126908468, -0.292156360634725, 0.292156360634725, -0.121015126908468},
		{0.0977197537924274, -0.255833636800846, 0.316227766016838, -0.255833636800846, 0.0977197537924273, 0.0977197537924281, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924271, 0.0977197537924277, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924270, 0.0977197537924279, -0.255833636800847, 0.316227766016838, -0.255833636800846, 0.0977197537924268},
		{0.0738219058991409, -0.205373505467450, 0.292156360634725, -0.315252941349890, 0.269628494389924, -0.165228553867128, 0.0248109445657181, 0.121015126908469, -0.240461479767574, 0.307490367669328, -0.307490367669328, 0.240461479767574, -0.121015126908468, -0.0248109445657177, 0.165228553867129, -0.269628494389925, 0.315252941349890, -0.292156360634724, 0.205373505467449, -0.0738219058991401},
		{0.0494689214077114, -0.143564401525505, 0.223606797749979, -0.281761002650515, 0.312334477467278, -0.312334477467278, 0.281761002650515, -0.223606797749979, 0.143564401525505, -0.0494689214077107, -0.0494689214077114, 0.143564401525505, -0.223606797749979, 0.281761002650515, -0.312334477467278, 0.312334477467278, -0.281761002650515, 0.223606797749978, -0.143564401525503, 0.0494689214077104},
		{0.0248109445657177, -0.0738219058991411, 0.121015126908468, -0.165228553867129, 0.205373505467450, -0.240461479767575, 0.269628494389924, -0.292156360634725, 0.307490367669328, -0.315252941349890, 0.315252941349890, -0.307490367669328, 0.292156360634725, -0.269628494389923, 0.240461479767575, -0.205373505467449, 0.165228553867128, -0.121015126908468, 0.0738219058991390, -0.0248109445657143}
	};

	std::vector<std::vector<double>> basisVectorCompare = dsp::fft::calcCosineBasisVectors<double>(nBasisVectors);
	ASSERT_EQ(basisVectorCompare.size(), nBasisVectors);
	for (int i = 0; i < nBasisVectors; ++i)
	{
		ASSERT_EQ(basisVectorCompare[i].size(), nBasisVectors);
		for (int k = 0; k < nBasisVectors; ++k)
		{
			EXPECT_TRUE(std::abs(basisVectors[i][k] - basisVectorCompare[i][k]) < epsi);
		}
	}

	std::vector<double> yDoubleCompare = dsp::fft::dct<double>(xDouble, nBasisVectors);
	ASSERT_EQ(xDouble.size(), nBasisVectors);
	for (int i = 0; i < nBasisVectors; ++i)
	{
		EXPECT_TRUE(std::abs(yDoubleCompare[i] - yDouble[i]) < epsi);
	}

	/*
	std::vector<int> y1 = dsp::fft::dct<int>(VInt1, nBasisVectors);
	EXPECT_EQ(yTestTarget1, y1);
	y1 = dsp::fft::dct<int>(VInt2, nBasisVectors);
	EXPECT_EQ(yTestTarget2, y1);
	y1 = dsp::fft::dct<int>(VInt3, nBasisVectors);
	EXPECT_EQ(yTestTarget3, y1);
	

	std::vector<float> y2 = dsp::fft::dct<float>(Vfloat1, nBasisVectors);
	EXPECT_EQ(yTestTarget1, y2);
	y2 = dsp::fft::dct<float>(Vfloat2, nBasisVectors);
	EXPECT_EQ(yTestTarget2, y2);
	y2 = dsp::fft::dct<float>(Vfloat3, nBasisVectors);
	EXPECT_EQ(yTestTarget3, y2);
	*/




}