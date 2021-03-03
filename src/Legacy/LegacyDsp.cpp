// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************
#define _USE_MATH_DEFINES
#include "LegacyDsp.h"
#include <cmath>

// Freq. of the musical note C0 (octave 0)
static const double FREQUENCY_C0 = 16.35159783;


// ****************************************************************************
// Berechnet x mod y **********************************************************

int modulo(int x, int y)
{
	if (y < 1) { y = 1; }

	if (x >= 0)
	{
		return (x % y);
	}

	// x < 0 ************************************************
	return (y - ((-x) % y));
}

// ****************************************************************************
// Gibt die Signalenergie für die Signalpunkte zwischen startPos und **********
// startPos + numSamples-1 zurück. ********************************************

double getSignalEnergy(const Signal& signal, int startPos, int numSamples)
{
	int i;
	double energy = 0.0;
	double value;

	if (numSamples < 0) { numSamples = 0; }
	int endPos = startPos + numSamples - 1;

	for (i = startPos; i <= endPos; i++)
	{
		value = signal.x[modulo(i, signal.N)];
		energy += value * value;
	}

	return energy;
}

// ****************************************************************************
// Gibt die Signalenergie für die Signalpunkte zwischen startPos und **********
// startPos + numSamples-1 zurück. ********************************************

double getSignalEnergy(const Signal16& signal, int startPos, int numSamples)
{
	int i;
	double energy = 0.0;
	double value;

	if (numSamples < 0) { numSamples = 0; }
	int endPos = startPos + numSamples - 1;

	for (i = startPos; i <= endPos; i++)
	{
		value = signal.x[modulo(i, signal.N)];
		energy += value * value;
	}

	return energy;
}

// ****************************************************************************
// Berechnet die mittlere Signalleistung für die Signalpunkte zwischen startPos 
// und startPos+numSamples-1. *************************************************

double getMeanSignalPower(Signal& signal, int startPos, int numSamples)
{
	double energy = getSignalEnergy(signal, startPos, numSamples);
	if (numSamples < 1) { numSamples = 1; }
	double power = energy / numSamples;

	return power;
}


// ****************************************************************************
// Wandelt Signale für Real- und Imaginäranteil in Signale für Amplitude und **
// Phase um. ******************************************************************

void rectangularToPolar(ComplexSignal& s, int length)
{
	s.setMinLength(length);
	double re, im;

	for (int i = 0; i < length; i++)
	{
		re = s.re[i];
		im = s.im[i];
		s.re[i] = sqrt(re * re + im * im);
		s.im[i] = atan2(im, re);
	}
}

// ****************************************************************************
// Wandelt Signale für Amplitude und Phase in Signale für Real- und ***********
// Imaginäranteil um. *********************************************************

void polarToRectangular(ComplexSignal& s, int length)
{
	s.setMinLength(length);
	double mag, phase;

	for (int i = 0; i < length; i++)
	{
		mag = s.re[i];
		phase = s.im[i];
		s.re[i] = mag * cos(phase);
		s.im[i] = mag * sin(phase);
	}
}

// ****************************************************************************
// Ergänzt die negativen Frequenzen in einem Spektrum *************************

void generateNegativeFrequencies(ComplexSignal* spectrum)
{
	if (spectrum == NULL) { return; }
	int i;
	int N = spectrum->N;

	for (i = N / 2 + 1; i < N; i++)
	{
		spectrum->re[i] = spectrum->re[N - i];
		spectrum->im[i] = -spectrum->im[N - i];
	}
}


// ****************************************************************************
// Berechnet die reelle DFT durch einfache Korrelation ************************

void realDFT(Signal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize)
{
	// mind. N/2+1 Werte für den Real- und Imaginäranteil allokieren
	timeSignal.setMinLength(length);
	freqSignal.setMinLength(length / 2 + 1);

	int i, k;
	double angle;
	double re, im;
	int l2 = length / 2;

	for (k = 0; k <= l2; k++)
	{
		re = 0.0;
		im = 0.0;

		for (i = 0; i < length; i++)
		{
			angle = (2.0 * M_PI * k * i) / (double)length;
			re += timeSignal.x[i] * cos(angle);
			im -= timeSignal.x[i] * sin(angle); // Imaginäranteile mit neg. Vorzeichen, wg. Konsistenz mit komplexer DFT
		}

		if (normalize)
		{
			im = -im;
			re /= (double)l2;
			im /= (double)l2;
			if ((k == 0) || (k == l2)) { re /= 2; }
		}

		freqSignal.re[k] = re;
		freqSignal.im[k] = im;
	}
}

// ****************************************************************************
// Berechnet die reelle inverse DFT aus freqSignal durch Korrelation. *********

void realIDFT(ComplexSignal& freqSignal, Signal& timeSignal, int length, bool normalize)
{
	// mind. N/2+1 Werte für den Real- und Imaginäranteil allokieren
	timeSignal.setMinLength(length);
	freqSignal.setMinLength(length / 2 + 1);

	int i, k;
	double re, im;
	double angle;

	// Alle x[i] müssen auf Null gesetzt werden *************
	for (i = 0; i < length; i++) { timeSignal.x[i] = 0.0; }

	for (k = 0; k <= length / 2; k++)
	{
		// Amplitude der Sinus- und Cosinusschwingung aus dem *
		// Real- und Imaginäranteil berechnen *****************

		if (normalize)
		{
			im = -freqSignal.im[k] / (double)(length / 2);   // Imaginäranteil wird negiert !
			re = freqSignal.re[k] / (double)(length / 2);
			if ((k == 0) || (k == length / 2)) { re /= 2.0; }
		}
		else
		{
			re = freqSignal.re[k];
			im = freqSignal.im[k];
		}

		// Sinus- und Cosinusschwingung der akt. Harmonischen add.
		for (i = 0; i < length; i++)
		{
			angle = (2.0 * M_PI * k * i) / (double)length;
			timeSignal.x[i] += re * cos(angle) + im * sin(angle);
		}
	}
}

// ****************************************************************************
// Es wird die komplexe schnelle FT auf dem Signal s der Länge ****************
// N = 2^lengthExponent berechnet. Das Ergebnis der Länge N steht am Ende *****
// wieder in s. ***************************************************************

void complexFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
	int i, j, k;

	// Signale wenn nötig verlängern ************************
	int N = 1 << lengthExponent;
	s.setMinLength(N);

	// Bitumsortierung vornehmen ****************************
	double tr, ti;
	int nm1 = N - 1;
	int nd2 = N / 2;

	j = nd2;
	for (i = 1; i <= N - 2; i++)
	{
		if (i < j)
		{
			tr = s.re[j];
			ti = s.im[j];
			s.re[j] = s.re[i];
			s.im[j] = s.im[i];
			s.re[i] = tr;
			s.im[i] = ti;
		}
		k = nd2;

		while (k <= j)
		{
			j -= k;
			k /= 2;
		}
		j += k;
	}   // next i


	// Ein Schleifendurchlauf für jede Stufe ****************
	int l, le, le2;
	double ur, ui;
	double sr, si;
	int jm1, ip;

	for (l = 1; l <= lengthExponent; l++)
	{
		le = 1 << l;
		le2 = le / 2;
		ur = 1.0;
		ui = 0.0;
		sr = cos(M_PI / (double)le2);
		si = -sin(M_PI / (double)le2);

		for (j = 1; j <= le2; j++)
		{
			jm1 = j - 1;
			for (i = jm1; i <= nm1; i += le)
			{
				ip = i + le2;
				tr = s.re[ip] * ur - s.im[ip] * ui;
				ti = s.re[ip] * ui + s.im[ip] * ur;
				s.re[ip] = s.re[i] - tr;
				s.im[ip] = s.im[i] - ti;
				s.re[i] += tr;
				s.im[i] += ti;
			}   // next i

			tr = ur;
			ur = tr * sr - ui * si;
			ui = tr * si + ui * sr;
		}   // next j
	}   // next l

	// Die Ergebniswerte normalisieren ? ******************************
	if (normalize)
	{
		for (i = 0; i < N; i++)
		{
			s.re[i] /= (double)N;
			s.im[i] /= (double)N;
		}
	}
}

// ****************************************************************************
// Es wird die inverse komplexe FFT auf dem komplexen Signal s berechnet. *****
// Die Signallänge beträgt 2^lengthExponent. **********************************

void complexIFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
	int i, k;

	// Signal wenn nötig verlängern *************************
	int N = 1 << lengthExponent;
	s.setMinLength(N);

	for (k = 0; k < N; k++) { s.im[k] = -s.im[k]; }

	complexFFT(s, lengthExponent, normalize);

	for (i = 0; i < N; i++) { s.im[i] = -s.im[i]; }
}

// ****************************************************************************
// Berechnet die inverse DFT für ein komplexes Signal durch Korrelation. ******

void complexIDFT(ComplexSignal& freqSignal, ComplexSignal& timeSignal, int length, bool normalize)
{
	int i, k;
	double angle;
	double S, C;

	freqSignal.setMinLength(length);

	// Das Frequenzbereichssignal auf Null setzen *********************
	timeSignal.reset(length);

	// Alle Zeitbereichsabtastwerte durchlaufen *********************

	for (i = 0; i < length; i++)
	{
		// Alle Frequenzen durchlaufen ************************************
		for (k = 0; k < length; k++)
		{
			angle = (2.0 * M_PI * k * i) / (double)length;
			S = sin(angle);
			C = cos(angle);
			timeSignal.re[i] += freqSignal.re[k] * C - freqSignal.im[k] * S;
			timeSignal.im[i] += freqSignal.im[k] * C + freqSignal.re[k] * S;
		}

		if (normalize)
		{
			timeSignal.re[i] /= length;
			timeSignal.im[i] /= length;
		}
	}
}


// ****************************************************************************
// Berechnet die DFT für ein komplexes Signal durch Korrelation. **************

void complexDFT(ComplexSignal& timeSignal, ComplexSignal& freqSignal, int length, bool normalize)
{
	int i, k;
	double angle;
	double S, C;

	timeSignal.setMinLength(length);

	// Das Frequenzbereichssignal auf Null setzen *********************
	freqSignal.reset(length);

	// Alle Frequenzen durchlaufen ************************************

	for (k = 0; k < length; k++)
	{
		// Alle Zeitbereichsabtastwerte durchlaufen *********************
		for (i = 0; i < length; i++)
		{
			angle = (2.0 * M_PI * k * i) / (double)length;
			S = sin(angle);
			C = cos(angle);
			freqSignal.re[k] += timeSignal.re[i] * C + timeSignal.im[i] * S;
			freqSignal.im[k] += timeSignal.im[i] * C - timeSignal.re[i] * S;
		}

		if (normalize)
		{
			freqSignal.re[k] /= length;
			freqSignal.im[k] /= length;
		}
	}
}


// ****************************************************************************
/// Returns the smallest exponent e for which windowLength <= 2^e.
/// For example, for windowLength = 500, e = 9, because 2^9 = 512.
// ****************************************************************************

int getFrameLengthExponent(int windowLength_pt)
{
	int e = 1;
	while (((int)1 << e) < windowLength_pt)
	{
		e++;
	}

	return e;
}


// ****************************************************************************
// Es wird die schnelle FT für ein reelles Signal berechnet, dass im Realteil
// von s steht. Der Imaginäranteil des Signals wird bei der Eingabe ignoriert.
// Das Ergebnis der Berechnung steht auch wieder in s (incl. der gespiegelten
// Frequenzen, wie bei der "normalen" FFT).
// N = 2^lengthExponent ist die Länge des Ein-/Ausgabesignals.
// Die Funktion ist ca. 30% schneller als die FFT mit komplexem Zeitsignal.
// ****************************************************************************

void realFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
	// Signale wenn nötig verlängern ************************
	int N = 1 << lengthExponent;
	s.setMinLength(N);

	int i, j, im, ip2, ipm, jm1, ip;

	// Die geraden und ungeraden Punkte trennen *************
	for (i = 0; i < N / 2; i++)
	{
		s.re[i] = s.re[2 * i];
		s.im[i] = s.re[2 * i + 1];
	}

	// Die normale FFT für N/2 Punkte berechnen *************
	complexFFT(s, lengthExponent - 1, false);     // Hier noch keine Normalisierung durchführen

	// Gerade/Ungerade Frequenzbereichszerlegung ************
	int nm1 = N - 1;
	int nd2 = N / 2;
	int n4 = (N / 4) - 1;

	for (i = 1; i <= n4; i++)
	{
		im = nd2 - i;
		ip2 = i + nd2;
		ipm = im + nd2;
		s.re[ip2] = (s.im[i] + s.im[im]) / 2;
		s.re[ipm] = s.re[ip2];
		s.im[ip2] = -(s.re[i] - s.re[im]) / 2;
		s.im[ipm] = -s.im[ip2];

		s.re[i] = (s.re[i] + s.re[im]) / 2;
		s.re[im] = s.re[i];
		s.im[i] = (s.im[i] - s.im[im]) / 2;
		s.im[im] = -s.im[i];
	}   // next i

	s.re[(N * 3) / 4] = s.im[N / 4];
	s.re[nd2] = s.im[0];
	s.im[(N * 3) / 4] = 0.0;
	s.im[nd2] = 0.0;
	s.im[N / 4] = 0.0;
	s.im[0] = 0.0;

	// Die letzte Stufe der FFT vervollständigen ************

	int l = lengthExponent;
	int le = 1 << l;
	int le2 = le / 2;
	double ur = 1.0;
	double ui = 0.0;
	double sr = cos(M_PI / (double)le2);
	double si = -sin(M_PI / (double)le2);
	double ti;
	double tr;

	for (j = 1; j <= le2; j++)
	{
		jm1 = j - 1;
		for (i = jm1; i <= nm1; i += le)
		{
			ip = i + le2;
			tr = s.re[ip] * ur - s.im[ip] * ui;
			ti = s.re[ip] * ui + s.im[ip] * ur;
			s.re[ip] = s.re[i] - tr;
			s.im[ip] = s.im[i] - ti;
			s.re[i] += tr;
			s.im[i] += ti;
		}   // next i

		tr = ur;
		ur = tr * sr - ui * si;
		ui = tr * si + ui * sr;
	}   // next j

	// Die Ergebniswerte normalisieren ? ******************************
	if (normalize)
	{
		for (i = 0; i < N; i++)
		{
			s.re[i] /= (double)N;
			s.im[i] /= (double)N;
		}
	}
}


// ****************************************************************************
// Es wird die schnelle inverse FT für reelle Signale berechnet.
// N = 2^lengthExponent ist die Länge der IDFT, und das Eingabesignal ist
// vom Index 0 bis zum Index N/2 mit den Spektralkoeffizienten gefüllt.
// Die restlichen Werte (N/2+1 .. N-1) werden symmetrisch ergänzt.
// Als Ausgabe enthält s.re[] das reelle Zeitsignal und s.im[] enthält Nullen.
// ****************************************************************************

void realIFFT(ComplexSignal& s, int lengthExponent, bool normalize)
{
	// Signale wenn nötig verlängern ************************
	int N = 1 << lengthExponent;
	s.setMinLength(N);

	int i, k;

	// Macht den Frequenzbereich symmetrisch ****************
	for (k = N / 2 + 1; k < N; k++)
	{
		s.re[k] = s.re[N - k];
		s.im[k] = -s.im[N - k];
	}

	// Real- und Imaginäranteil zusammenaddieren ************
	for (k = 0; k < N; k++) { s.re[k] += s.im[k]; }

	realFFT(s, lengthExponent, false);    // hier noch keine Normalisierung vornehmen

	// Postprocessing ***************************************
	for (i = 0; i < N; i++)
	{
		s.re[i] = s.re[i] + s.im[i];
		s.im[i] = 0.0;
	}

	// Die Ergebniswerte normalisieren ? ******************************
	if (normalize)
	{
		for (i = 0; i < N; i++) { s.re[i] /= (double)N; }
	}
}

// ****************************************************************************
// Berechnet ein Fenster der Länge l ******************************************

void getWindow(Signal& window, int length, WindowType type)
{
	int i;
	window.reset(length);

	// Ein einfaches Rechteckfenster ************************
	if (type == RECTANGULAR_WINDOW)
	{
		for (i = 0; i < length; i++) { window.x[i] = 1.0; }
	}

	// Ein Hamming-Window ***********************************
	else
		if (type == HAMMING_WINDOW)
		{
			for (i = 0; i < length; i++)
			{
				window.x[i] = 0.54 - 0.46 * cos((2.0 * M_PI * (double)i) / (double)(length - 1));
			}
		}

	// Die rechte Hälfte eines Hamming-Fensters *************
		else
			if (type == RIGHT_HALF_OF_HAMMING_WINDOW)
			{
				for (i = 0; i < length; i++)
				{
					window.x[i] = 0.54 - 0.46 * cos(M_PI + (M_PI * (double)i) / (double)(length - 1));
				}
			}

	// Die rechte Hälfte eines Hamming-Fensters *************
			else
				if (type == LEFT_HALF_OF_HAMMING_WINDOW)
				{
					for (i = 0; i < length; i++)
					{
						window.x[i] = 0.54 - 0.46 * cos((M_PI * (double)i) / (double)(length - 1));
					}
				}

	// Die rechte Hälfte eines Hann-Fensters ****************
				else
					if (type == RIGHT_HALF_OF_HANN_WINDOW)
					{
						for (i = 0; i < length; i++)
						{
							window.x[i] = 0.5 - 0.5 * cos(M_PI + (M_PI * (double)i) / (double)(length - 1));
						}
					}

	// Eine Gaußglocke **************************************
					else
						if (type == GAUSS_WINDOW)
						{
							const double yEdge = 0.01;
							double s = (double)(length * length) / (4.0 * log(yEdge));
							for (i = 0; i < length; i++)
							{
								window.x[i] = exp(((i - length / 2) * (i - length / 2)) / s);
							}
						}

	// Ein sonstiges Fenster ********************************
						else
						{
							for (i = 0; i < length; i++) { window.x[i] = 1.0; }
						}
}


// ----------------------------------------------------------------------------
// LPC-Funktionalität ---------------------------------------------------------
// ----------------------------------------------------------------------------


// ****************************************************************************
// Berechnet die LPC-Koeffizienten. coeff[0] wird auf 1 gesetzt. **************
// N gibt die Anz. Koeffizienten ausser coeff[0] an ! *************************

void getLPCCoefficients(const double* signal, int numSamples, double* coeff, int N)
{
	const int MAX_COEFF = 256;

	int i, j, p;
	double r[MAX_COEFF];
	double alpha[MAX_COEFF];
	double beta[MAX_COEFF];
	double z[MAX_COEFF];      // Reflektionskoeffizienten
	double E, q;

	// maximal 255 Koeffizienten erlauben *******************
	if (N > MAX_COEFF - 1) { N = MAX_COEFF - 1; }

	// Die lastCoeffIndex+1 Werte der Autokorellationsfolge nach r[]
	// berechnen ********************************************

	for (i = 0; i <= N; i++)
	{
		r[i] = 0.0;
		for (j = 0; j < numSamples - i; j++) { r[i] += signal[j] * signal[j + i]; }
	}

	// Der Levinson-Durbin-Algorithmus **********************

	E = r[0];
	alpha[0] = 1.0;
	z[0] = 0.0;

	for (p = 1; p <= N; p++)
	{
		q = 0.0;
		for (i = 0; i < p; i++) { q += alpha[i] * r[p - i]; }
		if (E == 0.0) { E = 0.0001; }
		z[p] = -q / E;
		alpha[p] = 0.0;
		for (i = 0; i <= p; i++) { beta[i] = alpha[i] + z[p] * alpha[p - i]; }
		for (i = 0; i <= p; i++) { alpha[i] = beta[i]; }

		E = E * (1.0 - z[p] * z[p]);
	}

	coeff[0] = 1;
	for (i = 1; i <= N; i++) { coeff[i] = -alpha[i]; }
}

// ****************************************************************************
// Berechnet das LPC-Restsignal nach residual[] durch inverse Filterung. ******

void getLPCResidual(const double* signal, double* residual, long l, const double* coeff, long N)
{
	long i, j;
	for (i = 0; i < l; i++)
	{
		residual[i] = signal[i];
		for (j = 1; j <= N; j++)
		{
			if (i - j >= 0) { residual[i] -= signal[i - j] * coeff[j]; }
		}
	}
}


// ****************************************************************************
// Berechnet das ursprüngliche Signal aus dem residual und den Koeffizienten. *

void predictSignal(double* signal, const double* residual, long l, const double* coeff, long N)
{
	long i, j;

	for (i = 0; i < l; i++)
	{
		signal[i] = residual[i];
		for (j = 1; j <= N; j++)
		{
			if (i - j >= 0) { signal[i] += coeff[j] * signal[i - j]; }
		}
	}
}

// ****************************************************************************
// Überführt die Koeffizienten des Prädiktorpolynoms 1-a1*(z^-1) - a2*(z^-2) - 
// ... -aN*(z^-N) zum Zwecke der Nullsetzung und durch Multiplikation mit z^N 
// in die Koeffizienten des Polynoms der Form a0*z^N + a1*z^(N-1) + ... + aN .
// Dazu müssen lediglich die Koeff. 1 .. N negiert werden. ********************

void LPCToPolynomCoefficients(double* LPCCoeff, double* polynomCoeff, long N)
{
	polynomCoeff[0] = LPCCoeff[0];
	long i;
	for (i = 1; i <= N; i++)
	{
		polynomCoeff[i] = -LPCCoeff[i];
	}
}

// ----------------------------------------------------------------------------
// Funktionen zur Nullstellenberechnung ---------------------------------------
// ----------------------------------------------------------------------------

// ****************************************************************************
// Gibt die (u.U. komplexen) Nullstellen der quadr. Gleichung x^2+beta*x+gamma 
// zurück. ********************************************************************

void getSquareRoots(double beta, double gamma, Complex& x0, Complex& x1)
{
	double re = -0.5 * beta;
	double root = 0.25 * beta * beta - gamma;

	// 2 komplexe Nullstellen *******************************
	if (root <= 0.0)
	{
		root = -root;
		double im = sqrt(root);
		x0 = Complex(re, im);
		x1 = Complex(re, -im);
	}
	// 2 reelle Nullstellen *********************************
	else
	{
		root = sqrt(root);
		x0 = Complex(re + root, 0.0);
		x1 = Complex(re - root, 0.0);
	}
}

// ****************************************************************************
// Gibt den Wert eines Polynoms an der Stelle x zurück. ***********************

Complex getPolynomValue(double* a, long N, Complex x)
{
	Complex result(0.0, 0.0);
	Complex power(1.0, 0.0);  //  Der Faktor x^n für hier x^0

	long i;

	for (i = N; i >= 0; i--)
	{
		result += a[i] * power;
		power *= x;
	}

	return result;
}


// ****************************************************************************
// Berechnet simulatan alle Quadratfaktoren des Polynoms N-ten Grades. Falls N 
// ungerade ist, wird es durch das Hinzufügen der neuen Nullstelle bei x=0 um *
// eins erhöht. ***************************************************************

void getPolynomRoots(double* a, int& N, Complex* roots)
{
	const long MAX_M = 128;   // max. Anzahl Quadratfaktoren

	int l, j, i, m;
	bool fertig = false;
	int E;                    // Zähler für noch ungenaue Quadratfaktoren
	double p, q;              // Koeffizienten des akt. Quadratvektors
	double c[2];              // Koeffizienten des Restglieds im Hornerschema
	double temp;
	double u, w;
	double S, T, Sm, Tm, SNew, TNew;
	double h, k;              // Die Korrekturen
	double D;                 // Die Determinante
	double startAngle[MAX_M]; // Wo liegt jeweils der 1. Näherungswert auf dem Einheitskreis ?


	// Falls der Polynomgrad N ungerade ist, dann eine weitere Nullstelle
	// durch Multiplikation des Polynoms mit x hinzufügen ***

	if ((N & 1) == 1)
	{
		N++;
		a[N] = 0.0;
	}

	int M = N / 2;       // Anzahl der Quadratfaktoren

	// Faktor[i] = x^2 + beta[i]*x + gamma[i]
	double beta[MAX_M];         // Hier werden die Indizes 1 .. M benutzt
	double gamma[MAX_M];

	// Die komplexen Einheitswurzeln als Startnäherung vorgeben
	for (m = 1; m <= M - 1; m++)
	{
		startAngle[m] = (M_PI * m) / (double)M;
		beta[m] = 2.0 * cos(startAngle[m]);
		gamma[m] = 1.0;
	}
	startAngle[M] = 0.0;
	beta[M] = 0.0;
	gamma[M] = -1;

	// relative Genauigkeitsschranke ************************
	double epsilon = 0.0001;
	double epsilon1 = 2.0 * epsilon;
	double epsilon2 = 2.0 * epsilon;
	double E1 = 0.0;
	double E2 = 0.0;
	const long MAX = 32;    // max. Anzahl der Iterationsschritte für alle Quadratfaktoren

	fertig = false;
	l = 1;

	while ((l <= MAX) && (!fertig))
	{
		E = 0;                // Zähler für die korrigierten Quadratfaktoren
		j = 1;

		// Alle M Quadratfaktoren durchlaufen
		while (j <= M)
		{
			p = beta[j];
			q = gamma[j];

			// Doppelzeiliges Hornerschema zur Bestimmung von c[0] und c[1]
			c[0] = a[0];
			c[1] = a[1] - p * a[0];
			for (i = 2; i <= N; i++)
			{
				temp = c[1];
				c[1] = a[i] - p * c[1] - q * c[0];
				c[0] = temp;
			}
			c[1] = c[1] + p * c[0];

			// Im ersten Durchlauf von l die Ungenauigkeiten aufsummieren
			if (l == 1) { E2 += fabs(c[0]) + fabs(c[1]); }

			// Weitere Korrektur ist noch nötig
			if (fabs(c[0]) + fabs(c[1]) >= epsilon2)
			{
				u = -0.5 * p;
				w = u * u - q;
				S = a[0];
				T = 0;

				// Berechnung des Produktes als S + v*T
				for (m = 1; m <= M; m++)
				{
					if (m != j)
					{
						Tm = beta[m] - p;
						Sm = u * Tm + gamma[m] - q;
						SNew = S * Sm + w * T * Tm;
						TNew = S * Tm + T * Sm;
						S = SNew;
						T = TNew;
					}
				}

				D = S * S - T * T * w;

				// Determinante ist nahe Null => Startwerte abändern und nochmal!
				if (fabs(D) < epsilon)
				{
					startAngle[j] += 0.012345;  // Winkel auf komplexem Zahlenkreis um ein paar Grad erhöhen
					beta[j] = 2.0 * cos(startAngle[j]);
					gamma[j] = 1.0;
					// j wird nicht erhöht !
				}
				else
					// Die Korrekturen für den j-ten Quadratfaktor können berechnet werden
				{
					h = (c[0] * (S - u * T) - T * c[1]) / D;
					k = (c[1] * (S + u * T) + c[0] * T * q) / D;
					beta[j] += h;
					gamma[j] += k;

					if (fabs(h) + fabs(k) >= epsilon1) { E++; }

					if (l == 1)
					{
						E1 += fabs(h) + fabs(k);
						epsilon1 = (epsilon * E1) / M;
						epsilon2 = (epsilon * E2) / M;
					}

					j++;    // Auf zum nächsten Quadratfaktor
				}
			}   // weitere Korrektur war nötig
			else { j++; }

		}   // Durchlauf der Quadratfaktoren

		if (E == 0) { fertig = true; }
		l++;

	}   // Iterationsschritte


	// ******************************************************
	// Aus den M beta[j] und gamma[j] die Nullstellen berechnen

	int nextRoot = 0;
	Complex comp[2];
	Complex re, wurzel;

	for (j = 1; j <= M; j++)
	{
		getSquareRoots(beta[j], gamma[j], comp[0], comp[1]);  // Nullstellen des Quadratfaktors
		roots[nextRoot++] = comp[0];
		roots[nextRoot++] = comp[1];
	}

	/*
	  {
		char st[512];
		sprintf(st, "Genauigkeit erreicht nach %d Schleifendurchläufen", l);
		MessageBox(NULL, st, "getPolynomRoots()", MB_OK);
	  }
	*/
}

// ****************************************************************************
// Berechnet simulatan alle Quadratfaktoren des Polynoms N-ten Grades. Falls N 
// ungerade ist, wird es durch das Hinzufügen der neuen Nullstelle bei x=0 um *
// eins erhöht. Zurückgegeben werden lediglich die reellen Nullstellen des ****
// Polynoms. ******************************************************************

void getRealPolynomRoots(double* a, int& N, double* roots, int& numRealRoots)
{
	const int MAX_M = 128;   // max. Anzahl Quadratfaktoren

	int l, j, i, m;
	bool fertig = false;
	int E;                    // Zähler für noch ungenaue Quadratfaktoren
	double p, q;              // Koeffizienten des akt. Quadratvektors
	double c[2];              // Koeffizienten des Restglieds im Hornerschema
	double temp;
	double u, w;
	double S, T, Sm, Tm, SNew, TNew;
	double h, k;              // Die Korrekturen
	double D;                 // Die Determinante
	double startAngle[MAX_M]; // Wo liegt jeweils der 1. Näherungswert auf dem Einheitskreis ?


	// Falls der Polynomgrad N ungerade ist, dann eine weitere Nullstelle
	// durch Multiplikation des Polynoms mit x hinzufügen ***

	if ((N & 1) == 1)
	{
		N++;
		a[N] = 0.0;
	}

	int M = N / 2;       // Anzahl der Quadratfaktoren

	// Faktor[i] = x^2 + beta[i]*x + gamma[i]
	double beta[MAX_M];         // Hier werden die Indizes 1 .. M benutzt
	double gamma[MAX_M];

	// Die komplexen Einheitswurzeln als Startnäherung vorgeben
	for (m = 1; m <= M - 1; m++)
	{
		startAngle[m] = (M_PI * m) / (double)M;
		beta[m] = 2.0 * cos(startAngle[m]);
		gamma[m] = 1.0;
	}
	startAngle[M] = 0.0;
	beta[M] = 0.0;
	gamma[M] = -1;

	// relative Genauigkeitsschranke ************************
	double epsilon = 0.0001;
	double epsilon1 = 2.0 * epsilon;
	double epsilon2 = 2.0 * epsilon;
	double E1 = 0.0;
	double E2 = 0.0;
	const long MAX = 32;    // max. Anzahl der Iterationsschritte für alle Quadratfaktoren

	fertig = false;
	l = 1;

	while ((l <= MAX) && (!fertig))
	{
		E = 0;                // Zähler für die korrigierten Quadratfaktoren
		j = 1;

		// Alle M Quadratfaktoren durchlaufen
		while (j <= M)
		{
			p = beta[j];
			q = gamma[j];

			// Doppelzeiliges Hornerschema zur Bestimmung von c[0] und c[1]
			c[0] = a[0];
			c[1] = a[1] - p * a[0];
			for (i = 2; i <= N; i++)
			{
				temp = c[1];
				c[1] = a[i] - p * c[1] - q * c[0];
				c[0] = temp;
			}
			c[1] = c[1] + p * c[0];

			// Im ersten Durchlauf von l die Ungenauigkeiten aufsummieren
			if (l == 1) { E2 += fabs(c[0]) + fabs(c[1]); }

			// Weitere Korrektur ist noch nötig
			if (fabs(c[0]) + fabs(c[1]) >= epsilon2)
			{
				u = -0.5 * p;
				w = u * u - q;
				S = a[0];
				T = 0;

				// Berechnung des Produktes als S + v*T
				for (m = 1; m <= M; m++)
				{
					if (m != j)
					{
						Tm = beta[m] - p;
						Sm = u * Tm + gamma[m] - q;
						SNew = S * Sm + w * T * Tm;
						TNew = S * Tm + T * Sm;
						S = SNew;
						T = TNew;
					}
				}

				D = S * S - T * T * w;

				// Determinante ist nahe Null => Startwerte abändern und nochmal!
				if (fabs(D) < epsilon)
				{
					startAngle[j] += 0.012345;  // Winkel auf komplexem Zahlenkreis um ein paar Grad erhöhen
					beta[j] = 2.0 * cos(startAngle[j]);
					gamma[j] = 1.0;
					// j wird nicht erhöht !
				}
				else
					// Die Korrekturen für den j-ten Quadratfaktor können berechnet werden
				{
					h = (c[0] * (S - u * T) - T * c[1]) / D;
					k = (c[1] * (S + u * T) + c[0] * T * q) / D;
					beta[j] += h;
					gamma[j] += k;

					if (fabs(h) + fabs(k) >= epsilon1) { E++; }

					if (l == 1)
					{
						E1 += fabs(h) + fabs(k);
						epsilon1 = (epsilon * E1) / M;
						epsilon2 = (epsilon * E2) / M;
					}

					j++;    // Auf zum nächsten Quadratfaktor
				}
			}   // weitere Korrektur war nötig
			else { j++; }

		}   // Durchlauf der Quadratfaktoren

		if (E == 0) { fertig = true; }
		l++;

	}   // Iterationsschritte


	// ******************************************************************
	// Aus den M beta[j] und gamma[j] die (reellen) Nullstellen berechnen

	numRealRoots = 0;
	double re, wurzel;

	for (j = 1; j <= M; j++)
	{
		re = -0.5 * beta[j];
		wurzel = 0.25 * beta[j] * beta[j] - gamma[j];

		if (wurzel >= 0.0)
		{
			wurzel = sqrt(wurzel);
			roots[numRealRoots++] = re + wurzel;
			roots[numRealRoots++] = re - wurzel;
		}
	}

	/*
	  {
		char st[512];
		sprintf(st, "Genauigkeit erreicht nach %d Schleifendurchläufen", l);
		MessageBox(NULL, st, "getPolynomRoots()", MB_OK);
	  }
	*/
}

// ****************************************************************************
// ****************************************************************************

double hertzToSemitones(double freq_Hz)
{
	if (freq_Hz < 1.0) { freq_Hz = 1.0; }
	return 12.0 * log(freq_Hz / FREQUENCY_C0) / log(2.0);
}

// ****************************************************************************
// ****************************************************************************

double semitonesToHertz(double freq_st)
{
	return FREQUENCY_C0 * pow(2, freq_st / 12.0);
}


// ****************************************************************************
