// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "LegacySignal.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

// ----------------------------------------------------------------------------
// Implementation of complex signals. 
// ----------------------------------------------------------------------------

// ****************************************************************************
// Constructor.
// ****************************************************************************

ComplexSignal::ComplexSignal(int length)
{
  re = NULL;
  im = NULL;
  N = 0;
  if (length > 0) { reset(length); }
}

// ****************************************************************************
// Destructor.
// ****************************************************************************

ComplexSignal::~ComplexSignal()
{
  dispose();
}

// ****************************************************************************
// Reset the signal to a new length an initialize with zeros.
// ****************************************************************************

void ComplexSignal::reset(int length)
{
  if (N != length)
  {
    if (re != NULL) { delete[] re; }
    if (im != NULL) { delete[] im; }
    N = length;
    re = NULL;
    im = NULL;
    if (N > 0) 
    { 
      re = new double[N]; 
      im = new double[N]; 
    }
  }

  if (N > 0) { setZero(); }
}

// ****************************************************************************
// Clear the signal.
// ****************************************************************************

void ComplexSignal::dispose()
{
  if (re != NULL)
  {
    delete[] re;
    re = NULL;
  }
  if (im != NULL)
  {
    delete[] im;
    im = NULL;
  }

  N = 0;
}

// ****************************************************************************
// Set all signal values to zero.
// ****************************************************************************

void ComplexSignal::setZero()
{
  int i;
  for (i=0; i < N; i++)
  {
    re[i] = im[i] = 0.0;
  }
}

// ****************************************************************************
// The length of the signal is changed to newLength. The content of the current
// signal is preserved (When newLength < current length, only the first 
// newLength values are preserved). New additional values are set to zero.
// ****************************************************************************

void ComplexSignal::setNewLength(int newLength)
{
  if (newLength != N)
  {
    //ComplexSignal copy = *this;
    ComplexSignal copy(N);
  	for (int i = 0; i < N; i++)
  	{
        copy.re[i] = this->re[i];
        copy.im[i] = this->im[i];
  	}
    copy.N = N;
    this->reset(newLength);
  
    int length = newLength;
    if (copy.N < length) { length = copy.N; }

    // Werte rüberkopieren ********************************
    memcpy(re, copy.re, length*sizeof(double));
    memcpy(im, copy.im, length*sizeof(double));
  }
}

// ****************************************************************************
// When the current signal length < minLength, then it is set to minLength and
// the current signal values are copied. After calling this function the signal
// length is at least newLength.
// ****************************************************************************

void ComplexSignal::setMinLength(int minLength)
{
  if (N < minLength) { this->setNewLength(minLength); }
}

// ****************************************************************************
// Make sure that (0 <= index < N).
// ****************************************************************************

void ComplexSignal::limitIndex(int& index)
{
  if (N > 0)
  {
    if (index < 0)
      { index = N - ((-index) % N); }
    else
      { index = index % N; }
  }
}


// ****************************************************************************
// Set a value.
// ****************************************************************************

void ComplexSignal::setValue(int pos, Complex value)
{
  if (N > 0)
  {
    limitIndex(pos);
    re[pos] = value.real();
    im[pos] = value.imag();
  }
}

// ****************************************************************************
// Set a value.
// ****************************************************************************

void ComplexSignal::setValue(int pos, double newRe, double newIm)
{
  if (N > 0)
  {
    limitIndex(pos);
    re[pos] = newRe;
    im[pos] = newIm;
  }
}

// ****************************************************************************
// Return different values at a specific position of the signal.
// ****************************************************************************

Complex ComplexSignal::getValue(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return Complex(re[pos], im[pos]);
  }
  return Complex(0.0, 0.0);
}

// ****************************************************************************

double ComplexSignal::getMagnitude(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return sqrt(re[pos]*re[pos] + im[pos]*im[pos]);
  }
  return 0.0;
}

// ****************************************************************************

double ComplexSignal::getPhase(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return atan2(im[pos], re[pos]);
  }
  return 0.0;
}

// ****************************************************************************

double ComplexSignal::getRealPart(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return re[pos];
  }
  return 0.0;
}

// ****************************************************************************

double ComplexSignal::getImaginaryPart(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return im[pos];
  }
  return 0.0;
}

// ****************************************************************************
// Operators.
// ****************************************************************************

void ComplexSignal::operator=(ComplexSignal& s)
{
  reset(s.N);
  if (re != NULL) 
  { 
    memcpy(re, s.re, N*sizeof(double)); 
  }
  if (im != NULL) 
  { 
    memcpy(im, s.im, N*sizeof(double)); 
  }
}

// ****************************************************************************

void ComplexSignal::operator+=(ComplexSignal& s)
{
  setMinLength(s.N);

  // Werte addieren ***************************************
  for (int i=0; i < s.N; i++) 
  { 
    re[i]+= s.re[i]; 
    im[i]+= s.im[i]; 
  }
}

// ****************************************************************************

void ComplexSignal::operator*=(ComplexSignal& s)
{
  double newRe, newIm;
  setMinLength(s.N);

  // komplexe Zahlen multiplizieren ***********************
  for (int i=0; i < s.N; i++) 
  { 
    newRe = re[i]*s.re[i] - im[i]*s.im[i];
    newIm = re[i]*s.im[i] + im[i]*s.re[i];
    re[i] = newRe;
    im[i] = newIm; 
  }
}

// ****************************************************************************

void ComplexSignal::operator*=(double factor)
{
  // mit konstantem Faktor multiplizieren *****************
  for (int i=0; i < N; i++) 
  { 
    re[i]*= factor;
    im[i]*= factor; 
  }
}

// ****************************************************************************
