// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __LEGACY__SIGNAL_H__
#define __LEGACY__SIGNAL_H__

#include <complex>
#include <memory.h>

typedef std::complex<double> Complex;

template<class ElementType> class TemplateSignal;

typedef TemplateSignal<double>       Signal;
typedef TemplateSignal<signed short> Signal16;
typedef TemplateSignal<int>          Signal32;


// ****************************************************************************
// Template for signals with different element types.
// ****************************************************************************

template<class ElementType> class TemplateSignal
{
  public:
    int N;               // Length of the signal
    ElementType* x;      // The values

    TemplateSignal(int length = 0);
    ~TemplateSignal();

    void reset(int length);
    void dispose();
    void setZero();
    void setNewLength(int newLength);
    void setMinLength(int minLength);
    void limitIndex(int& index);
    void setValue(int pos, ElementType value);
    ElementType getValue(int pos);
    void getMinMax(ElementType &min, ElementType &max);

    void writeTo(TemplateSignal<ElementType>& s, int startPos, bool wrap=true);

    // Operators **************************************************************

    void operator=(TemplateSignal<ElementType>& s);
    void operator+=(TemplateSignal<ElementType>& s);
    void operator*=(TemplateSignal<ElementType>& s);
    void operator*=(double factor);
};


// ****************************************************************************
// Class for signals with complex values.
// ****************************************************************************

class ComplexSignal
{
  public:
    int N;          // Length of the signal 
    double* re;     // Real value parts
    double* im;     // Imaginary value parts

    ComplexSignal(int length = 0);
    ~ComplexSignal();

    void reset(int length);
    void dispose();
    void setZero();
    void setNewLength(int newLength);
    void setMinLength(int minLength);
    void limitIndex(int& index);
    void setValue(int pos, Complex value);
    void setValue(int pos, double newRe, double newIm);
    
    Complex getValue(int pos);
    double getMagnitude(int pos);
    double getPhase(int pos);
    double getRealPart(int pos);
    double getImaginaryPart(int pos);

    // Operators **************************************************************

    void operator=(ComplexSignal& s);
    void operator+=(ComplexSignal& s);
    void operator*=(ComplexSignal& s);
    void operator*=(double factor);
};


// ----------------------------------------------------------------------------
// Template-Signal Implementation. 
// Must be in the same file as the declaration!
// ----------------------------------------------------------------------------

// ****************************************************************************
// Constructor.
// ****************************************************************************

template<class ElementType> TemplateSignal<ElementType>::TemplateSignal(int length)
{
  x = NULL;
  N = 0;
  if (length > 0) { reset(length); }
}

// ****************************************************************************
// Destructor.
// ****************************************************************************

template<class ElementType> TemplateSignal<ElementType>::~TemplateSignal()
{
  dispose();
}

// ****************************************************************************
// Initialization with zeros.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::reset(int length)
{
  if (N != length)
  {
    if (x != NULL) { delete[] x; }
    N = length;
    x = NULL;
    if (N > 0) { x = new ElementType[N]; }
  }

  // Mit Nullen initialisieren ****************************
  if (N > 0) { setZero(); }
}    

// ****************************************************************************
// Free the memory.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::dispose()
{
  if (x != NULL)
  {
    delete[] x;
    x = NULL;
  }
  N = 0;
}

// ****************************************************************************
// All signal values to zeros.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::setZero()
{
  int i;
  for (i=0; i < N; i++) { x[i] = (ElementType)0; }
}

// ****************************************************************************
// The length of the signal is changed to newLength. The content of the current
// signal is preserved (When newLength < current length, only the first 
// newLength values are preserved). New additional values are set to zero.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::setNewLength(int newLength)
{
  if (newLength != N)
  {
    //TemplateSignal<ElementType> copy = *this;
    TemplateSignal<ElementType> copy(N);
  	for (int i = 0; i < N; i++)
  	{
        copy.x[i] = this->x[i];
  	}
    copy.N = N;
    this->reset(newLength);
  
    int length = newLength;
    if (copy.N < length) { length = copy.N; }

    // Werte r�berkopieren ********************************
    memcpy(x, copy.x, length*sizeof(ElementType));
  }
}

// ****************************************************************************
// When the current signal length < minLength, then it is set to minLength and
// the current signal values are copied. After calling this function the signal
// length is at least newLength.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::setMinLength(int minLength)
{
  if (N < minLength) { this->setNewLength(minLength); }
}

// ****************************************************************************
// Make sure that (0 <= index < N).
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::limitIndex(int& index)
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
// Returns the signal value at the position pos. A modulo calculation is 
// applied, when pos is out of the signal range.
// ****************************************************************************

template<class ElementType> ElementType TemplateSignal<ElementType>::getValue(int pos)
{
  if (N > 0)
  {
    limitIndex(pos);
    return x[pos];
  }
  
  return 0;
}

// ****************************************************************************
// Returns the minimum and maximum values in the signal.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::getMinMax(ElementType &min, ElementType &max)
{
  if (N < 1)
  {
    min = 0;
    max = 0;
  }

  min = max = x[0];
  int i;
  for (i=1; i < N; i++)
  {
    if (x[i] > max) { max = x[i]; }
    if (x[i] < min) { min = x[i]; }
  }
}

// ****************************************************************************
// Sets the signal value at the position pos. A modulo calculation is 
// applied, when pos is out of the signal range.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::setValue(int pos, ElementType value)
{
  if (N > 0)
  {
    limitIndex(pos);
    x[pos] = value;
  }
}

// ****************************************************************************
// Writes the entire signal to the position startPos in the signal s. The 
// corresponding points in s are overwritten. wrap tells, wheather points
// bejond the limit of s are wrapped to the beginning of s.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::writeTo(TemplateSignal<ElementType>& s, int startPos, bool wrap)
{
  if (s.N > 0)
  {
    s.limitIndex(startPos);

    for (int i=0; i < N; i++)
    {
      s.x[startPos++] = x[i];
      if (startPos >= s.N) 
      { 
        if (wrap) 
          { startPos = 0; }
        else
          { i = N; }
      }
    }
  }

}

// ****************************************************************************
// This signal becomes the copy of the signal on the right side of the "=".
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::operator=(TemplateSignal<ElementType>& s)
{
  reset(s.N);
  if (x != NULL) 
  { 
    memcpy(x, s.x, N*sizeof(ElementType)); 
  }
}

// ****************************************************************************
// This signal is added value by value with s. The resulting signal length is
// at least the length of signal s.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::operator+=(TemplateSignal<ElementType>& s)
{
  setMinLength(s.N);
  for (int i=0; i < s.N; i++) { x[i]+= s.x[i]; }
}

// ****************************************************************************
// This signal is multiplied value by value with s. The resulting signal length 
// is at least the length of signal s.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::operator*=(TemplateSignal<ElementType>& s)
{
  setMinLength(s.N);
  for (int i=0; i < s.N; i++) { x[i]*= s.x[i]; }
}

// ****************************************************************************
// All signal values are multiplied by factor.
// ****************************************************************************

template<class ElementType> void TemplateSignal<ElementType>::operator*=(double factor)
{
	for (int i=0; i < N; i++) { x[i] = (ElementType)(x[i]*factor); }
}

// ****************************************************************************

#endif
