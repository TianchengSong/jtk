/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song of CGG.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/LocalCorrelationFilter.java

The original documentation for Local correlation filter.

/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/**
 * Local cross-correlation of two arrays with seamless overlapping windows.
 * Given two input arrays f and g and a specified lag, this filter computes 
 * an output array c of local cross-correlation coefficients, one for each
 * sample in the input arrays f and g.
 * <p>
 * Two types of cross-correlation are implemented. Both types can be 
 * normalized to obtain cross-correlation coefficients with magnitudes 
 * that do not exceed one. The normalization varies, depending on the 
 * type of cross-correlation.
 * <p>
 * <em>Simple</em> cross-correlation computes an array of products 
 * h[j] = f[j]*g[j+lag] and then filters this array of products with a 
 * window. The resulting correlation cfg[k,lag] is not symmetric with 
 * respect to lag; cfg[k,-lag] = cgf[k-lag,lag] != cgf[k,lag]. For
 * simple cross-correlation, normalization scale factors vary with lag
 * and should be applied before picking correlation peaks.
 * <p>
 * <em>Symmetric</em> cross-correlation computes an array of products
 * h[j] = f[j-lag/2]*g[j+lag/2] and therefore requires interpolation
 * between samples for odd lags. (For efficiency, we interpolate the 
 * products h, not the inputs f and g.) The resulting correlation is 
 * symmetric with respect to lag; cfg[k,lag] = cgf[k,-lag]. Moreover,
 * when inputs f and g are the same, each local auto-correlation has a
 * Fourier transform (a power spectrum) that is positive-semidefinite.
 * For symmetric cross-correlation, normalization scale factors do not 
 * vary with lag, and therefore need not be applied before picking 
 * correlation peaks.
 * <p>
 * Two correlation windows are implemented: Gaussian and rectangular.
 * Gaussian windows should be used for most applications. Rectangular
 * windows are provided primarily for comparison, because they are so
 * often used by others.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.08.11
 */

/******************************************************************************
 *  HAMPSON - RUSSELL SOFTWARE SERVICES LTD.    CALGARY, ALBERTA, CANADA
 ******************************************************************************/
/*!
 *  \class    LocalCorrelationFilter
 * 
 *  \brief    Dave Hale's local correlation filter.
 *            The naming of the variable does not follow the HRS convention.
 *            The reason is to follow Dave Hale's coding as much as possible so
 *            any future algorithmic change can be made easily.
 *            In the first implementation, only the single dimension methods
 *            are implemented.
 *            All pointer operations are changed to std::vector.
 *
 ******************************************************************************/
/* CREATOR:     Tiancheng Song                     DATE: Nov, 2014
  *****************************************************************************
 * REVISIONS:
 * AUTHOR:                                         DATE: 
 * DESCRIPTION:  
  *****************************************************************************
*/
#ifndef HrsLocalCorrelationFilter_include
#define HrsLocalCorrelationFilter_include

#include <string>
#include <vector>
#include "HrsRecursiveGaussianFilter.h"
using namespace std;
class GaussianFilter;

class LocalCorrelationFilter {
public:
  /**
   * Cross-correlations types.
   */
  enum Type {
    SIMPLE,
    SYMMETRIC
  };

  /**
   * Cross-correlations windows.
   */
   enum Window {
    GAUSSIAN,
    RECTANGLE
  };
  LocalCorrelationFilter(Type type, Window window, double sigma);
  LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2);
  LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2, double sigma3);
  void setInputs(const vector<float>& f, const vector<float>& g);
  void setInputs(const vector<vector<float> >& f, const vector<vector<float> >& g);
  void setInputs(const vector<vector<vector<float> > >& f, const vector<vector<vector<float> > >& g);
  void correlate(int lag, vector<float>& c);
  void correlate(int lag1, int lag2, vector<vector<float> >& c);
  void correlate(int lag1, int lag2, int lag3, vector<vector<vector<float> > >& c);
  void normalize(int lag, vector<float>& c) ;
  void normalize(int lag1, int lag2, vector<vector<float> >& c) ;
  void normalize(int lag1, int lag2, int lag3, vector<vector<vector<float> > >& c);
  vector<float> unbias(const vector<float>& f) ;
  vector<vector<float> > unbias(vector<vector<float> >& f);
  vector<vector<vector<float> > > unbias(const vector<vector<vector<float> > >& f);
  ///////////////////////////////////////////////////////////////////////////
private:
  static const float S1;
  static const float S2;
  static const float S3;
  static const float S4;
  static const float S[8];
  Window _window; // window for correlations
  Type _type;     // correlation type
  double _sigma1,_sigma2,_sigma3; // window half-widths
  int _dimension; // dimension of input arrays; 0 if no inputs
  int _n1,_n2,_n3; // array lengths
  vector<vector<vector<float> > >  _f,_g; // inputs f and g; by reference
  vector<vector<vector<vector<float> > > > _s; // normalization scale factors

  void correlate(int lag, const vector<float>& f, const vector<float>& g, vector<float>& c);
  void correlate(
    int lag1, int lag2, const vector<vector<float> >& f, const vector<vector<float> >& g, vector<vector<float> >& c);
  void correlate(
    int lag1, int lag2, int lag3, const vector<vector<vector<float> > >& f, const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& c); 
  void updateNormalize();
  static void shift(const vector<float>& f, vector<float>& g);
  static void shift1(const vector<vector<float> >& f, vector<vector<float> >& g);
  static void shift2(const vector<vector<float> >& f, vector<vector<float> >& g);
  static void shift1(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g);

  static void shift2(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g);
  static void shift3(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g);
  void checkDimension(int dimension);
  void checkDimensions(vector<float>& c);
  void checkDimensions(vector<vector<float> >& c);
  void checkDimensions(vector<vector<vector<float> > >& c) ;

  static void sub(const vector<float>&a,const vector<float>&b, vector<float>&c);	
  static void sub(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c);
  static void sub(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c);
  static void sqrtArray(const vector<float>&a, vector<float>&c);
  static void sqrtArray(const vector<vector<float> >& a, vector<vector<float> >& c);
  static void sqrtArray(const vector<vector<vector<float> > >& a, vector<vector<vector<float> > >& c);
  static void div(const vector<float>&a,const vector<float>&b, vector<float>&c);	
  static void div(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c);
  static void div(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c);
  static void div(const float a,const vector<float>&b, vector<float>&c);	
  static void div(const float a,const vector<vector<float> >& b, vector<vector<float> >& c);
  static void div(const float a,const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c);
  static void multiplyArray(const vector<float>&a,const vector<float>&b, vector<float>&c);	
  static void multiplyArray(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c);
  static void multiplyArray(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c);
  // This interface makes it easier to implement different windows.
  GaussianFilter *_f1, *_f2, *_f3; // filters used to implement windows
};

class GaussianFilter{
  public:
    GaussianFilter(double sigma) {
      _rgf = new HrsRecursiveGaussianFilter(sigma);
    }
    void apply(const vector<float>& x, vector<float>& y) {
      _rgf->apply0(x,y);
    }
    void apply1(const vector<vector<float> >& x, vector<vector<float> >& y) {
      _rgf->apply0X(x,y);
    }
    void apply2(const vector<vector<float> >& x, vector<vector<float> >& y) {
      _rgf->applyX0(x,y);
    }
    void apply1(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
      _rgf->apply0XX(x,y);
    }
    void apply2(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
      _rgf->applyX0X(x,y);
    }
    void apply3(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
      _rgf->applyXX0(x,y);
    }
  private:
    HrsRecursiveGaussianFilter* _rgf;
  };

#endif
