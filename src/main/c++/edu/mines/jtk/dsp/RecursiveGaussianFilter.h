/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/RecursiveGaussianFilter.java

The original documentation for Recursive Gaussian Filter.

/****************************************************************************
Copyright (c) 2005, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/**
 * Recursive implementation of a Gaussian filter and derivatives. Filters 
 * include the 0th, 1st, and 2nd derivatives. The impulse response of the 
 * 0th-derivative smoothing filter is infinitely long, and is approximately 
 * h[n] = 1.0/(sqrt(2*PI)*sigma)*exp(-0.5*(n*n)/(sigma*sigma)). Here,
 * sigma denotes the standard width of the Gaussian.
 * <p>
 * For large filter widths sigma, this recursive implementation can be
 * much more efficient than convolution with a truncated Gaussian.
 * Specifically, if the Gaussian is truncated for |n| &gt; 4*sigma, then
 * this recursive implementation requires 2/sigma of the multiplications 
 * required by convolution. In other words, for sigma &gt; 2, this
 * recursive implementation should be more efficient than convolution.
 * <p>
 * For any application of this filter, input and output arrays may be the 
 * same array. When the filter cannot be applied in-place, intermediate
 * arrays are constructed internally.
 * <p>
 * This filter implements two different methods for approximating 
 * with difference equations a Gaussian filter and its derivatives.
 * <p>
 * The first method is that of Deriche, R., 1993, Recursively implementing 
 * the Gaussian and its derivatives: INRIA Research Report, number 1893. 
 * Deriche's method is used for small widths sigma, for which it is most 
 * accurate. 
 * <p>
 * The second method is that of van Vliet, L.J., Young, I.T., and Verbeek, 
 * P.W., 1998, Recursive Gaussian derivative filters, Proceedings of the 
 * 14th International Conference on Pattern Recognition, IEEE Computer 
 * Society Press. The parallel implementation used here yields zero-phase 
 * impulse responses without the end effects caused by the serial (cascade) 
 * poles-only implementation recommended by van Vliet, et al. This 
 * second method is used for large widths sigma.
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.02.12
*/

#ifndef RecursiveGaussianFilter_include
#define RecursiveGaussianFilter_include

#include "CDouble.h"
#include "Recursive2ndOrderFilter.h"
#include <string>
#include <vector>
#include <stdint.h> // for uint64_t
#include <sstream>  // for std::ostringstream

using namespace std;

class Filter;
class Cdouble;
class VanVlietFilter;
class DericheFilter;
class Recursive2ndOrderFilter;

class Filter {
public:
  virtual void applyN(int nd, const vector<float>& x, vector<float>& y){};
  virtual void applyXN(int nd, const vector< vector<float> >& x, vector< vector<float> >& y){};
  void applyNX(int nd, const vector< vector<float> >& x, vector< vector<float> >& y);
  void applyNXX(int nd, const vector<vector<vector<float> > >& x,  vector<vector<vector<float> > >& y);
  void applyXNX(int nd, const vector<vector<vector<float> > >& x,  vector<vector<vector<float> > >& y);
  void applyXXN(int nd, const vector<vector<vector<float> > >& x,  vector<vector<vector<float> > >& y);
};
/******************************* INCLUDE FILES ********************************/
class  RecursiveGaussianFilter // __HRS_DECLSPEC_SEISUTIL
{
public:
  // The method used to design the Gaussian filter.  
  enum Method 
  { DERICHE,    
    VAN_VLIET  
  };

  RecursiveGaussianFilter(double sigma, Method method);
  RecursiveGaussianFilter(double sigma);
  ~RecursiveGaussianFilter();
  void apply0(const vector<float>& x,  vector<float>& y);
  void apply1(const vector<float>& x,  vector<float>& y);
  void apply2(const vector<float>& x,  vector<float>& y);
  void apply0X(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply1X(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply2X(const vector< vector<float> >& x, vector< vector<float> >& y);
  void applyX0(const vector< vector<float> >& x, vector< vector<float> >& y);
  void applyX1(const vector< vector<float> >& x, vector< vector<float> >& y);
  void applyX2(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply0XX(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply1XX(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply2XX(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyX0X(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyX1X(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyX2X(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyXX0(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyXX1(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void applyXX2(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply00(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply10(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply01(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply11(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply20(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply02(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply000(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply100(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply010(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply001(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply110(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply101(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply011(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply200(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply020(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply002(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  static void checkArrays(const vector<float>& x, const vector<float>& y);
  static void checkArrays(const vector< vector<float> >& x, const vector< vector<float> >& y);
  static void checkArrays(const vector<vector<vector<float> > >& x, const vector<vector<vector<float> > >& y);
  static bool sameArrays(const vector<float>& x, const vector<float>& y);
  static bool sameArrays(const vector< vector<float> >& x, const vector< vector<float> >& y);
  float getSigma(){return mySigma;};
  void setSigma(float sigma){mySigma = sigma;};
private:
  Filter* _filter;
  float mySigma;
};

class __HRS_DECLSPEC_SEISUTIL DericheFilter: public Filter {
public:
    DericheFilter(double sigma);
    virtual void applyN(int nd, const vector<float>& x, vector<float>& y);
    virtual void applyXN(int nd, const vector< vector<float> >& x, vector< vector<float> >& y);

    // Coefficients computed using Deriche's method. These coefficients
    // were computed for sigma = 100 and 0 <= x <= 10*sigma = 1000,
    // using the Mathematica function FindFit. The coefficients have
    // roughly 10 digits of precision.
    // 0th derivative.
private:
  static double a00; 
  static double a10;
  static double b00;
  static double b10;
  static double c00;
  static double c10;
  static double w00;
  static double w10;
  // 1st derivative.
  static double a01;
  static double a11;
  static double b01;
  static double b11;
  static double c01; 
  static double c11;
  static double w01;
  static double w11;
  // 2nd derivative.
  static double a02;
  static double a12;
  static double b02;
  static double b12;
  static double c02;
  static double c12;
  static double w02;
  static double w12;
  //
  static double a0[];
  static double a1[];
  static double b0[];
  static double b1[];
  static double c0[];
  static double c1[];
  static double w0[];
  static double w1[];

  vector<float> _n0,_n1,_n2,_n3; // numerator coefficients
  vector<float> _d1,_d2,_d3,_d4; // denominator coefficients

  /**
    * Makes Deriche's numerator and denominator coefficients.
    */
  void makeND(double sigma);
  /**
    * Scales numerator filter coefficients to normalize the filters.
    * For example, the sum of the 0th-derivative filter coefficients
    * should be 1.0. The scale factors are computed from finite-length
    * approximations to the impulse responses of the three filters.
    */
  void scaleN(double sigma);
};

///////////////////////////////////////////////////////////////////////////
class VanVlietFilter: public Filter{
public:
  VanVlietFilter(double sigma);
  virtual void applyN(int nd, const vector<float>&x,  vector<float>&y) ;
  virtual void applyXN(int nd, const vector< vector<float> >& x, vector< vector<float> >& y);

private:
  // Poles (inverses) for 4th-order filters published by van Vliet, et al.
  vector< vector<Cdouble> > POLES;
  vector<vector<vector<Recursive2ndOrderFilter> > > _g;

  void makeG(double sigma);
  Recursive2ndOrderFilter makeFilter(double b0, double b1, double b2, double a1, double a2);
  Cdouble gr(int nd, Cdouble polej, const vector<Cdouble>& poles, double gain) ;
  static vector<Cdouble> adjustPoles(double sigma, const vector<Cdouble>& poles) ;
  static double computeGain(vector<Cdouble>&  poles);
  static double computeSigma(double sigma, const vector<Cdouble>& poles);
};

#endif
