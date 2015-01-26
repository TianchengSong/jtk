/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song of CGG.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/LocalShiftFinder.java

The original documentation for LocalShiftFinder.
*/
/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/**
 * Estimates displacement vector fields for two images. For example, given 
 * two 2-D images f(x1,x2) and g(x1,x2), a shift finder estimates two vector 
 * components of displacement u1(x1,x2) and u2(x1,x2) such that 
 * f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * <p>
 * Like the images f and g, the components of displacement are sampled
 * functions of coordinates x1 and x2. That is, displacements may vary 
 * from sample to sample. The components u1 and u2 represent displacements 
 * in the x1 and x2 coordinate directions, respectively.
 * <p>
 * This shift finder estimates each component of displacement using local
 * cross-correlations. For each image sample, the estimated shift is that 
 * which yields the maximum correlation coefficient. This coefficient is
 * found by quadratic interpolation of correlation functions sampled at
 * integer lags.
 * <p>
 * The peak (maximum) correlation coefficient may be used to measure 
 * quality of an estimated shift. However, because a correlation function 
 * may have more than one peak (local maxima), a better measure of quality 
 * may be the difference between the coefficients for the correlation peak 
 * and next highest peak. Both the peak coefficient and this difference may 
 * be computed with the shifts.
 * <p>
 * Methods are provided to find and compensate for each component of shift 
 * sequentially. As each component is found, that component can be removed 
 * from the image g before estimating another component. For example, again 
 * for 2-D images f(x1,x2) and g(x1,x2), we might first estimate u1(x1,x2). 
 * If we then compute an image h(x1,x2) = g(x1+u1(x1,x2),x2), we can use
 * f(x1,x2) and h(x1,x2) to estimate u2(x1,x2). By repeating this process
 * sequentially, we obtain estimates for both u1(x1,x2) and u2(x1,x2) such
 * that f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * <p>
 * Methods are also provided to whiten 2-D and 3-D images before estimating
 * displacement vectors. This (spectral) whitening improves estimates of
 * displacements parallel to image features that may be otherwise poorly
 * resolved. Whitening is performed with local prediction error filters
 * computed from local auto-correlations.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.11.18
 */

#ifndef HrsLocalShiftFinder_include
#define HrsLocalShiftFinder_include

#include "HrsDeclSpec.h"
#include "LocalCorrelationFilter.h"
#include <string>
#include <vector>
#include "SincInterpolator.h"

using namespace std;

class LocalCorrelationFilter;
class HrsRecursiveGaussianFilter;
class SincInterpolator;

class __HRS_DECLSPEC_SEISUTIL LocalShiftFinder {

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has the
   * same correlation window half-width for all dimensions.
   * @param sigma the correlation window half-width; must not be less than 1.
   */
public:
  LocalShiftFinder(double sigma);
  LocalShiftFinder(double sigma1, double sigma2);
  LocalShiftFinder(double sigma1, double sigma2, double sigma3);
  void setInterpolateDisplacements(bool enable);
  void find1(int min1, int max1, const vector<float>& f, const vector<float>& g, vector<float>& u);
  void find1(int min1, int max1, const vector<float>& f, const vector<float>& g, 
    vector<float>& u, vector<float>& c, vector<float>& d);
  void find1(int min1, int max1, const vector<vector<float> >& f, 
    const vector<vector<float> >& g, vector<vector<float> >& u);
  void find2(int min2, int max2, const vector<vector<float> >& f, 
    const vector<vector<float> >& g, vector<vector<float> >& u) ;
  void find1(int min1, int max1, const vector<vector<vector<float> > >& f, 
    const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& u) ;
  void find2(int min2, int max2, const vector<vector<vector<float> > >& f, 
    const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& u);
  void find3(int min3, int max3, const vector<vector<vector<float> > >& f, 
    const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& u);
  void shift1(const vector<float>& du, vector<float>& u1, vector<float>& h) ;
  void shift1(const vector<vector<float> >& du, vector<vector<float> >& u1, vector<vector<float> >& u2, vector<vector<float> >& h) ;
  void shift2(const vector<vector<float> >& du, vector<vector<float> >& u1, vector<vector<float> >& u2, vector<vector<float> >& h);
  void shift1(const vector<vector<vector<float> > >& du, vector<vector<vector<float> > >& u1, 
    vector<vector<vector<float> > >& u2, vector<vector<vector<float> > >& u3,vector<vector<vector<float> > >& h); 
  void shift2(const vector<vector<vector<float> > >& du, vector<vector<vector<float> > >& u1, vector<vector<vector<float> > >& u2, 
    vector<vector<vector<float> > >& u3,vector<vector<vector<float> > >& h);
  void shift3(const vector<vector<vector<float> > >& du, vector<vector<vector<float> > >& u1, vector<vector<vector<float> > >& u2, 
    vector<vector<vector<float> > >& u3,vector<vector<vector<float> > >& h);
  void whiten(const vector<vector<float> >& f,  vector<vector<float> >& g);
  void whiten(double sigma, const vector<vector<float> >& f, vector<vector<float> >& g);
  void whiten(const vector<vector<vector<float> > >& f,  vector<vector<vector<float> > >& g);
  void whiten(double sigma, const vector<vector<vector<float> > >& f,  vector<vector<vector<float> > >& g);

private:
  LocalCorrelationFilter* _lcfSimple;
  SincInterpolator* _si;
  bool _interpolateDisplacements;

  void findShifts(
    int min, int max, const vector<float>& f, const vector<float>& g, vector<float>& u, vector<float>& c, vector<float>& d);
  void findShifts(
    int min, int max, const vector<float>& f, const vector<float>& g, vector<float>& u);
  void findShifts(
    int dim, int min, int max, const vector<vector<float> >& f, const vector<vector<float> >& g, vector<vector<float> >& u);
  void findShifts(
    int dim, int min, int max, const vector<vector<vector<float> > >& f, const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& u);
};

#endif
