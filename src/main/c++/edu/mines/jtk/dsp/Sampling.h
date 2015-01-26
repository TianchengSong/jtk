/*
The algorithm contained in the .cxx and .h are used with permission from Dave Hale 
of Colorado School of Mines and is under the terms of Common Public License v1.0.

The code is a straight port from java to C++ by Tiancheng Song of CGG.
The complete java source code can be downloaded from Dave Hale's website:
https://github.com/dhale/jtk/blob/master/src/main/java/edu/mines/jtk/dsp/Sampling.java
*/

/****************************************************************************
Copyright (c) 2005, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

#ifndef Sampling_include
#define Sampling_include

#include <vector>
#include <cmath>

using namespace std;

class Sampling {
public:
  static double DEFAULT_TOLERANCE;
  Sampling(int n);
  Sampling(int n, double d, double f) ;
  Sampling(int n, double d, double f, double t);
  Sampling(vector<double>& v);
  Sampling(vector<double>& v, double t);
  int getCount() ;
  double getDelta();
  double getFirst();
  double getLast() ;
  double getValue(int i) ;
  vector<double> getValues();
  bool isUniform() ;
  bool isEquivalentTo(Sampling s);
  bool isCompatible(Sampling s) ;
  int indexOf(double x) ;
  int indexOfNearest(double x);
  double valueOfNearest(double x) ;
  bool isInBounds(int i);
  bool isInBounds(double x);
  bool isInBoundsExtended(double x) ;
  double getValueExtended(int i);
  int indexOfNearestExtended(double x);
  double valueOfNearestExtended(double x) ;
  int indexOfFloorExtended(double x);
  double normalizedDifferenceExtended(double x, int i);
  vector<int> overlapWith(Sampling s);
  Sampling mergeWith(Sampling s);
  Sampling shift(double s) ;
  Sampling prepend(int m);
  Sampling append(int m);
  Sampling decimate(int m);
  Sampling interpolate(int m);
private:
  int _n; // number of samples
  double _d; // sampling interval
  double _f; // value of first sample
  vector<double> _v; // array[n] of sample values; null, if uniform
  double _t; // sampling tolerance, as a fraction of _d
  double _td; // sampling tolerance _t multiplied by _d
  bool almostEqua;

  void init(int n, double d, double f, double t);
  double value(int i);
  bool almostEqual(double v1, double v2, double tiny) ;
  double tinyWith(Sampling s);
  int round(double r) {return (r > 0.0) ? (int)(r + 0.5) : (int)(r - 0.5);}
};
#endif
