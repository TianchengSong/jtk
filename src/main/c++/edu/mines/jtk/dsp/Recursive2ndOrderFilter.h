/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/HrsRecursive2ndOrderFilter.java

The original documentation for Recursive 2nd Order Filter.

/****************************************************************************
Copyright (c) 2005, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/**
 * Recursive 2nd-order filter. This filter solves a linear, 2nd-order, 
 * constant-coefficient difference equation in either forward or reverse
 * directions along any dimension of a 1-D, 2-D, or 3-D array.
 * <p>
 * Application of the filter in the forward direction computes
 * <pre>
 * y[i] = b0*x[i]+b1*x[i-1]+b2*x[i-2]-a1*y[i-1]-a2*y[i-2],
 * </pre>
 * for i = 0, 1, 2, ..., n-1, where x[i] = y[i] = 0 for i&lt;0.
 * Application of the filter in the reverse direction computes
 * <pre>
 * y[i] = b0*x[i]+b1*x[i+1]+b2*x[i+2]-a1*y[i+1]-a2*y[i+2],
 * </pre>
 * for i = n-1, n-2, ..., 0, where x[i] = y[i] = 0 for i&gt;=n.
 * @author Dave Hale, Colorado School of Mines
 * @version 2005.11.22
 */

#ifndef Recursive2ndOderFilter_include
#define Recursive2ndOderFilter_include

#include "CDouble.h"
#include <string>
#include <vector>

using namespace std;

class HrsRecursive2ndOrderFilter {
public:
  HrsRecursive2ndOrderFilter(
    float b0, float b1, float b2, float a1, float a2);
  HrsRecursive2ndOrderFilter(double pole, double zero, double gain);
  HrsRecursive2ndOrderFilter(
    Cdouble pole1, Cdouble pole2,
    Cdouble zero1, Cdouble zero2,
    double gain);
  HrsRecursive2ndOrderFilter();
  void applyForward(const vector<float>& x, vector<float>& y);
  void applyReverse(const vector<float>& x, vector<float>& y);
  void accumulateForward(const vector<float>& x, vector<float>& y);
  void accumulateReverse(const vector<float>& x, vector<float>& y);
  void apply1Forward(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply1Reverse(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply2Forward(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply2Reverse(const vector< vector<float> >& x, vector< vector<float> >& y);
  void accumulate1Forward(const vector< vector<float> >& x, vector< vector<float> >& y);
  void accumulate1Reverse(const vector< vector<float> >& x, vector< vector<float> >& y);
  void accumulate2Forward(const vector< vector<float> >& x, vector< vector<float> >& y);
  void accumulate2Reverse(const vector< vector<float> >& x, vector< vector<float> >& y);
  void apply1Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply1Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply2Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply2Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply3Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void apply3Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate1Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate1Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate2Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate2Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate3Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
  void accumulate3Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y);
private:
  float _b0,_b1,_b2,_a1,_a2; // filter coefficients
  static void checkArrays(const vector<float>& x, const vector<float>& y);
  static void checkArrays(const vector< vector<float> >& x, const vector< vector<float> >& y);
  static void checkArrays(const vector<vector<vector<float> > >& x, const vector<vector<vector<float> > >& y);
  void get2(int i2, const vector<vector<vector<float> > >& x, vector< vector<float> >& x2);
  void set2(int i2, const vector< vector<float> >& x2, vector<vector<vector<float> > >& x);
  void acc2(int i2, const vector< vector<float> >& x2, vector<vector<vector<float> > >& x);
};

#endif
