#include "gtest/gtest.h"

#include <numeric>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>


#include "dsp.h"
#include "dct_ref_data.h"

struct DspTest : ::testing::Test
{

};


TEST_F(DspTest, SignalMethods)
{
	auto s = dsp::Signal<double>(1000);
	EXPECT_TRUE(s.getSamplingRate_Hz() == 1000);
	s.getSamplingRate_Hz() = 2000;
	EXPECT_TRUE(s.getSamplingRate_Hz() == 2000);
	auto sr = &s.getSamplingRate_Hz();
	*sr = 3000;
	EXPECT_TRUE(s.getSamplingRate_Hz() == 3000);
}
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
	double epsi = 1e-5;

	std::vector<std::vector<long double>> basisVectorLongDoubleCompare = calcCosineBasisVectors<long double>(20);
	std::vector<std::vector<double>> basisVectorDoubleCompare = calcCosineBasisVectors<double>(20);
	std::vector<std::vector<float>> basisVectorFloatCompare = calcCosineBasisVectors<float>(20);

	// Test the basis vector calculation.
	ASSERT_EQ(basisVectorLongDoubleCompare.size(), 20);
	ASSERT_EQ(basisVectorDoubleCompare.size(), 20);
	ASSERT_EQ(basisVectorFloatCompare.size(), 20);
	for (int i = 0; i < 20; ++i)
	{
		ASSERT_EQ(basisVectorLongDoubleCompare[i].size(), 20);
		ASSERT_EQ(basisVectorDoubleCompare[i].size(), 20);
		ASSERT_EQ(basisVectorFloatCompare[i].size(), 20);
		for (int k = 0; k < 20; ++k)
		{
			EXPECT_NEAR(basisVectorsLongDouble[i][k], basisVectorLongDoubleCompare[i][k], epsi);
			EXPECT_NEAR(basisVectorsDouble[i][k], basisVectorDoubleCompare[i][k], epsi);
			EXPECT_NEAR(basisVectorsFloat[i][k], basisVectorFloatCompare[i][k], epsi);
		}
	}

	// Test the dct function when signal length == n.
	int n = 10;
	std::vector<long double> yLongDoubleCompare = dsp::fft::dct<long double>(xLongDouble, n, dsp::fft::dctType::dct2);
	std::vector<double> yDoubleCompare = dsp::fft::dct<double>(xDouble, n, dsp::fft::dctType::dct2);
	std::vector<float> yFloatCompare = dsp::fft::dct<float>(xFloat, n, dsp::fft::dctType::dct2);

	ASSERT_EQ(yLongDoubleCompare.size(), n);
	ASSERT_EQ(yDoubleCompare.size(), n);
	ASSERT_EQ(yFloatCompare.size(), n);
	for (int i = 0; i < n; ++i)
	{
		EXPECT_NEAR(yLongDoubleCompare[i], yLongDouble[i], epsi);
		EXPECT_NEAR(yDoubleCompare[i], yDouble[i], epsi);
		EXPECT_NEAR(yFloatCompare[i], yFloat[i], epsi);
	}

	// Test the dct function when signal length < n.
	n = 15;
	yLongDoubleCompare = dsp::fft::dct<long double>(xLongDouble, n, dsp::fft::dctType::dct2);
	yDoubleCompare = dsp::fft::dct<double>(xDouble, n, dsp::fft::dctType::dct2);
	yFloatCompare = dsp::fft::dct<float>(xFloat, n, dsp::fft::dctType::dct2);

	ASSERT_EQ(yLongDoubleCompare.size(), n);
	ASSERT_EQ(yDoubleCompare.size(), n);
	ASSERT_EQ(yFloatCompare.size(), n);
	for (int i = 0; i < n; ++i)
	{
		EXPECT_NEAR(yLongDoubleCompare[i], yLongDouble15[i], epsi);
		EXPECT_NEAR(yDoubleCompare[i], yDouble15[i], epsi);
		EXPECT_NEAR(yFloatCompare[i], yFloat15[i], epsi);
	}

	// Test the dct function when signal length > n.
	n = 5;
	yLongDoubleCompare = dsp::fft::dct<long double>(xLongDouble, n, dsp::fft::dctType::dct2);
	yDoubleCompare = dsp::fft::dct<double>(xDouble, n, dsp::fft::dctType::dct2);
	yFloatCompare = dsp::fft::dct<float>(xFloat, n, dsp::fft::dctType::dct2);

	ASSERT_EQ(yLongDoubleCompare.size(), n);
	ASSERT_EQ(yDoubleCompare.size(), n);
	ASSERT_EQ(yFloatCompare.size(), n);
	for (int i = 0; i < n; ++i)
	{
		EXPECT_NEAR(yLongDoubleCompare[i], yLongDouble5[i], epsi);
		EXPECT_NEAR(yDoubleCompare[i], yDouble5[i], epsi);
		EXPECT_NEAR(yFloatCompare[i], yFloat5[i], epsi);
	}
}