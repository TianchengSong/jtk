/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song of CGG.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/SincInterpolator.java

The original documentation for SincInterpolator.
*/
/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/**
 * A sinc interpolator for bandlimited uniformly-sampled functions y(x). 
 * Interpolators can be designed for any two of three parameters: maximum 
 * error (emax), maximum frequency (fmax) and maximum length (lmax). The 
 * parameter not specified is computed when an interpolator is designed.
 * <p>
 * Below the specified (or computed) maximum frequency fmax, the maximum 
 * interpolation error should be less than the specified (or computed) 
 * maximum error emax. For frequencies above fmax, interpolation error 
 * may be much greater. Therefore, sequences to be interpolated should 
 * be bandlimited to frequencies less than fmax.
 * <p>
 * The maximum length lmax of an interpolator is an even positive integer. 
 * It is the number of uniform samples required to interpolate a single 
 * value y(x). Ideally, the weights applied to each uniform sample are 
 * values of a sinc function. Although the ideal sinc function yields zero 
 * interpolation error for all frequencies up to the Nyquist frequency 
 * (0.5 cycles/sample), it has infinite length.
 * <p>
 * With recursive filtering, infinite-length approximations to the sinc
 * function are feasible and, in some applications, most efficient. When
 * the number of interpolated values is large relative to the number of 
 * uniform samples, the cost of recursive filtering is amortized over those 
 * many interpolated values, and can be negligible. However, this cost 
 * becomes significant when only a few values are interpolated for each 
 * sequence of uniform samples.
 * <p>
 * This interpolator is based on a <em>finite-length</em> approximation 
 * to the sinc function. The efficiency of finite-length interpolators 
 * like this one does not depend on the number of samples interpolated.
 * Also, this interpolator is robust in the presence of noise spikes, 
 * which affect only nearby samples.
 * <p>
 * Finite-length interpolators present a tradeoff between cost and accuracy.
 * Interpolators with small maximum lengths are most efficient, and those 
 * with high maximum frequencies and small maximum errors are most accurate.
 * <p>
 * When interpolating multiple values of y(x) from a single sequence of
 * uniformly sampled values, efficiency may be improved by using one of the
 * methods that enables specification of multiple x values at which to
 * interpolate.
 *
 * @author Dave Hale, Colorado School of Mines
 * @author Bill Harlan, Landmark Graphics
 * @version 2012.12.21
 */
#ifndef HrsSincInterpolator_include
#define HrsSincInterpolator_include

#include <vector>
#include <stdint.h>
#include <cmath>
#include "Sampling.h"
#include <assert.h>
#include <map>

#ifndef PI
#define PI (3.141592653589793)
#endif

using namespace std;

class KaiserWindow;
static inline uint64_t doubleToLongBits(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof bits);
    return bits;
}

class SincInterpolator {
 public:
   enum Extrapolation {
    ZERO, 
    CONSTANT
  };

  static SincInterpolator fromErrorAndLength(
    double emax, int lmax);
 
  static SincInterpolator fromErrorAndFrequency(
    double emax, double fmax);

  static SincInterpolator fromFrequencyAndLength(
    double fmax, int lmax);

  SincInterpolator() ;

  double getMaximumError();

  double getMaximumFrequency() ;

  int getMaximumLength();

  long getTableBytes();

 
  Extrapolation getExtrapolation();

  void setExtrapolation(Extrapolation extrap);

  float interpolate(
    int nxu, double dxu, double fxu, vector<float>& yu, double xi);
  
  void interpolate(
    int nxu, double dxu, double fxu, vector<float>& yu, 
    int nxi, vector<float>& xi, vector<float>& yi);

  void interpolate(
    int nxu, double dxu, double fxu, vector<float>& yu, 
    int nxi, double dxi, double fxi, vector<float>& yi);

  float interpolate(
    int nx1u, double dx1u, double fx1u, 
    int nx2u, double dx2u, double fx2u, 
    vector<vector<float> >& yu, double x1i, double x2i);


  float interpolate(
    int nx1u, double dx1u, double fx1u, 
    int nx2u, double dx2u, double fx2u, 
    int nx3u, double dx3u, double fx3u, 
    vector<vector<vector<float> > >& yu, double x1i, double x2i, double x3i);

  float interpolate(Sampling sxu, vector<float>& yu, double xi);
  void interpolate(
    Sampling sxu, vector<float>& yu, 
    Sampling sxi, vector<float>& yi);

  float interpolate(
    Sampling sx1u, Sampling sx2u,
    vector<vector<float> >& yu, double x1i, double x2i);

  
  float interpolate(
    Sampling sx1u, Sampling sx2u, Sampling sx3u,
    vector<vector<vector<float> > >& yu, double x1i, double x2i, double x3i);

  
  void interpolateComplex(
    int nxu, double dxu, double fxu, vector<float>& yu, 
    int nxi, double dxi, double fxi, vector<float>& yi);

  void interpolateComplex(
    int nxu, double dxu, double fxu, vector<float>& yu, 
    int nxi, vector<float>& xi, vector<float>& yi);

  void interpolateComplex(
    Sampling sxu, vector<float>& yu, 
    Sampling sxi, vector<float>& yi);

 
  void accumulate(
    double xa, float ya,
    int nxu, double dxu, double fxu, vector<float>& yu);
  
  void accumulate(
    int nxa, vector<float>& xa, vector<float>& ya,
    int nxu, double dxu, double fxu, vector<float>& yu);

 
  vector<vector<float> > getTable();

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Fraction of error due to windowing; remainder is due to table lookup.
 private:
   static double EWIN_FRAC;

  // Maximum table size, when maximum error and frequency are specified.
  static int NTAB_MAX;

  // Extrapolation method.
  Extrapolation _extrap;// = ZERO;

  // Table of sinc interpolation coefficients.

  int _lsinc; // length of sinc approximations
  int _nsinc; // number of sinc approximations
  double _dsinc; // sampling interval in table
  vector<vector<float> > _asinc; // array[nsinc][lsinc] of sinc approximations
  double _nsincm1; // nsinc-1
  int _ishift; // -lsinc-lsinc/2+1
  SincInterpolator(double emax,double fmax,int lmax);

  void init(double emax,double fmax,int lmax);

public:
  class Design {
  public:
    double emax;
    double fmax; 
    int lmax;
    Design(double emax, double fmax, int lmax) {
      this->emax = emax;
      this->fmax = fmax;
      this->lmax = lmax;
    }
    Design() {
      this->emax = 0.0;
      this->fmax = 0.0;
      this->lmax = 0;
    }
    int hashCode() {
      uint64_t lemax = doubleToLongBits(emax);
      uint64_t lfmax = doubleToLongBits(fmax);
      return (int)(lemax^(lemax >> 32)) ^
             (int)(lfmax^(lfmax >> 32)) ^
             lmax;
    }
    bool equals(Design& object) {      
      return this->emax==object.emax &&
             this->fmax==object.fmax &&
             this->lmax==object.lmax;
    }
  };

  /**
   * Table of sinc interpolator coefficients.
   */
   //class Table {
 struct Table {
   Design design; // here, all three design parameters are non-zero
   int lsinc,nsinc,nsincm1,ishift;
   double dsinc;
   vector<vector<float> > asinc;
 };
  static Table makeTable(Design& design); 
  static Table makeTable(int nsinc, int lsinc, KaiserWindow& kwin);
  static double sinc(double x);
  //static HashMap<Design,Table> _tables = new HashMap<Design,Table>();
  static map<Design,Table> _tables;// = new HashMap<Design,Table>();
  //map<Design,Table> _tables;// = new HashMap<Design,Table>();
  static Table getTable(double emax, double fmax, int lmax);
  Table _table; // with all fields cached below
  float interpolate(
    double xscale, double xshift, int nxum, int nxu, 
    vector<float>& yu, double x);

  void shift(
    int nxu, double dxu, double fxu, vector<float>& yu,
    int nxi,             double fxi, vector<float>& yi);
 
  void accumulate(
    double xscale, double xshift, int nxum,
    double x, float y, int nxu, vector<float>& yu);

  float interpolate(
    double x1scale, double x1shift, int nx1um, int nx1u,
    double x2scale, double x2shift, int nx2um, int nx2u,
    vector<vector<float> >& yu, double x1, double x2);

  float interpolate(
    double x1scale, double x1shift, int nx1um, int nx1u,
    double x2scale, double x2shift, int nx2um, int nx2u,
    double x3scale, double x3shift, int nx3um, int nx3u,
    vector<vector<vector<float> > >& yu, double x1, double x2, double x3);

  void interpolateComplex(
    double xscale, double xshift, int nxum, int nxu, 
    vector<float>& yu, int ix, double x, vector<float>& y);
};


class KaiserWindow {

  /**
   * Returns a Kaiser window with specified error and transition width.
   * @param error the maximum absolute error.
   * @param width the transition width.
   * @return the window.
   */
public:
  static KaiserWindow* fromErrorAndWidth(double error, double width) {
    assert(error>0.0);
    assert(error<1.0);
    assert(width>0.0);
    double a = -20.0*log10(error);
    double d = (a>21.0)?(a-7.95)/14.36:0.9222;
    double length = d/width;
    KaiserWindow* ptr = new KaiserWindow(error,width,length);
    return ptr;
  }

  /**
   * Returns a Kaiser window with specified error and window length.
   * @param error the maximum absolute error.
   * @param length the two-sided window length.
   * @return the window.
   */
   static KaiserWindow* fromErrorAndLength(double error, double length) {
    assert(error>0.0);
    assert(error<1.0);
    assert(length>0);
    double a = -20.0*log10(error);
    double d = (a>21.0)?(a-7.95)/14.36:0.9222;
    double width = d/length;
    KaiserWindow* ptr = new KaiserWindow(error,width,length);
    return ptr;
  }

  /**
   * Returns a Kaiser window with specified transition width and window length.
   * The product width*length cannot be less than one.
   * @param width the transition width
   * @param length the two-sided window length.
   * @return the window.
   */
   static KaiserWindow* fromWidthAndLength(double width, double length) {
    assert(width>0.0);
    assert(length>0);
    assert(width*length>=1.0);
    double d = width*length;
    double a = 14.36*d+7.95;
    double error = pow(10.0,-a/20.0);
    KaiserWindow* ptr = new KaiserWindow(error,width,length);
    return ptr;
  }

  /**
   * Returns the value of this Kaiser window function w(x) for specified x.
   * @param x the argument for which to evaluate w(x).
   * @return the value w(x).
   */
  double evaluate(double x) {
    double xx = x*x;
    return (xx<=_xxmax)?_scale*ino(_alpha*sqrt(1.0-xx/_xxmax)):0.0;
  }

  /**
   * Gets the maximum absolute error.
   * @return the maximum absolute error.
   */
  double getError() {
    return _error;
  }

  /**
   * Gets the two-sided window length.
   * @return the window length.
   */
  double getLength() {
    return _length;
  }

  /**
   * Gets the transition width.
   * @return the transition width.
   */
  double getWidth() {
    return _width;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

private:
  double _error;
  double _width;
  double _length;
  double _alpha;
  double _scale;
  double _xxmax;

  KaiserWindow(double error, double width, double length) {
    _error = error;
    _width = width;
    _length = length;
    double a = -20.0*log10(_error);
    if (a<=21.0) {
      _alpha = 0.0;
    } else if (a<=50.0) {
      _alpha = 0.5842*pow(a-21.0,0.4)+0.07886*(a-21.0);
    } else {
      _alpha = 0.1102*(a-8.7);
    }
    _scale = 1.0/ino(_alpha);
    _xxmax = 0.25*_length*_length;
  }

  double ino(double x) {
    double s = 1.0;
    double ds = 1.0;
    double d = 0.0;
    do {
      d += 2.0;
      ds *= (x*x)/(d*d);
      s += ds;
    } while (ds>s*DBL_EPSILON);
    return s;
  }
};
#endif
