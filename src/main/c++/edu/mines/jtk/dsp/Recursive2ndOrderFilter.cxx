/*! \file Recursive2ndOrderFilter.cxx */

/*
 * C++ port from Java by Tiancheng Song
 * Please refer to the header file for complete documentation and licensing information.
 *
 * Member Functions:
 *       Recursive2ndOrderFilter::checkArrays
 *       Recursive2ndOrderFilter::get2
 *       Recursive2ndOrderFilter::set2
 *       Recursive2ndOrderFilter::acc2
 *       Recursive2ndOrderFilter::Recursive2ndOrderFilter
 *       Recursive2ndOrderFilter::applyForward
 *       Recursive2ndOrderFilter::accumulateForward
 *       Recursive2ndOrderFilter::accumulateReverse
 *       Recursive2ndOrderFilter::apply1Forward
 *       Recursive2ndOrderFilter::apply1Reverse
 *       Recursive2ndOrderFilter::apply2Forward
 *       Recursive2ndOrderFilter::apply2Reverse
 *       Recursive2ndOrderFilter::accumulate1Forward
 *       Recursive2ndOrderFilter::accumulate1Reverse
 *       Recursive2ndOrderFilter::accumulate2Forward
 *       Recursive2ndOrderFilter::accumulate2Reverse
 *       Recursive2ndOrderFilter::apply3Forward
 *       Recursive2ndOrderFilter::apply3Reverse
 *
 ******************************************************************************
*/
#include "Recursive2ndOrderFilter.h"
#include <assert.h>

/**
  * Constructs a recursive 2nd-order filter with specified coefficients.
  * If some of the coefficients are zero, the filter may be of only 1st
  * or even 0th order.
  * @param b0 a filter coefficient.
  * @param b1 a filter coefficient.
  * @param b2 a filter coefficient.
  * @param a1 a filter coefficient.
  * @param a2 a filter coefficient.
  */
Recursive2ndOrderFilter::Recursive2ndOrderFilter(
    float b0, float b1, float b2, float a1, float a2) 
{
  _b0 = b0;
  _b1 = b1;
  _b2 = b2;
  _a1 = a1;
  _a2 = a2;
}

/**
  * Constructs a recursive 2nd-order filter from pole, zero, and gain.
  * This filter is actually a 1st-order filter, because it has only
  * one (real) pole and zero.
  * @param pole the pole.
  * @param zero the zero.
  * @param gain the filter gain.
  */
Recursive2ndOrderFilter::Recursive2ndOrderFilter(double pole, double zero, double gain) {
  _b0 = (float)(gain);
  _b1 = (float)(-gain*zero);
  _a1 = (float)(-pole);
}
Recursive2ndOrderFilter::Recursive2ndOrderFilter()
{
}
/**
  * Constructs a recursive 2nd-order filter from poles, zeros, and gain.
  * The poles must be real or conjugate pairs; likewise for the zeros.
  * @param pole1 the 1st pole.
  * @param pole2 the2nd pole.
  * @param zero1 the 1st zero.
  * @param zero2 the 2nd zero.
  * @param gain the filter gain.
  */
Recursive2ndOrderFilter::Recursive2ndOrderFilter(
  Cdouble pole1, Cdouble pole2,
  Cdouble zero1, Cdouble zero2,
  double gain)
{
  assert(pole1.i==0.0 &&  pole2.i==0.0 ||
          pole2.r==pole1.r && -pole2.i==pole1.i);
  assert(zero1.i==0.0  &&  zero2.i==0.0 ||
          zero2.r==zero1.r && -zero2.i==zero1.i);
  _b0 = (float)(gain);
  _b1 = (float)(-(zero1.r+zero2.r)*gain);
  _b2 = (float)((zero1.times(zero2)).r*gain);
  _a1 = (float)(-(pole1.r+pole2.r));
  _a2 = (float)((pole1.times(pole2)).r);
}

/**
  * Applies this filter in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must have equal
  * lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::applyForward(const vector<float>& x, vector<float>& y) {
  checkArrays(x,y);
  int n = y.size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    float yim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi-_a1*yim1;
      y[i] = yi;
      yim1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    float yim1 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1-_a1*yim1;
      y[i] = yi;
      yim1 = yi;
      xim1 = xi;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1-_a1*yim1-_a2*yim2;
      y[i] = yi;
      yim2 = yim1;
      yim1 = yi;
      xim1 = xi;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim2 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b1*xim1+_b2*xim2-_a1*yim1-_a2*yim2;
      y[i] = yi;
      yim2 = yim1;
      yim1 = yi;
      xim2 = xim1;
      xim1 = xi;
    }
  }

  // General case.
  else { 
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim2 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1+_b2*xim2-_a1*yim1-_a2*yim2;
      y[i] = yi;
      yim2 = yim1;
      yim1 = yi;
      xim2 = xim1;
      xim1 = xi;
    }
  }
}

/**
  * Applies this filter in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must have equal
  * lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::applyReverse(const vector<float>& x, vector<float>& y) {
  checkArrays(x,y);
  int n = y.size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    float yip1 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi-_a1*yip1;
      y[i] = yi;
      yip1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    float xip1 = 0.0f;
    float yip1 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1-_a1*yip1;
      y[i] = yi;
      yip1 = yi;
      xip1 = xi;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    float xip1 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1-_a1*yip1-_a2*yip2;
      y[i] = yi;
      yip2 = yip1;
      yip1 = yi;
      xip1 = xi;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    float xip1 = 0.0f;
    float xip2 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b1*xip1+_b2*xip2-_a1*yip1-_a2*yip2;
      y[i] = yi;
      yip2 = yip1;
      yip1 = yi;
      xip2 = xip1;
      xip1 = xi;
    }
  }

  // General case.
  else { 
    float xip1 = 0.0f;
    float xip2 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1+_b2*xip2-_a1*yip1-_a2*yip2;
      y[i] = yi;
      yip2 = yip1;
      yip1 = yi;
      xip2 = xip1;
      xip1 = xi;
    }
  }
}

/**
  * Applies this filter in the forward direction, accumulating the output. 
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must have equal
  * lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulateForward(const vector<float>& x, vector<float>& y) {
  checkArrays(x,y);
  int n = y.size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    float yim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi-_a1*yim1;
      y[i] += yi;
      yim1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    float yim1 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1-_a1*yim1;
      y[i] += yi;
      yim1 = yi;
      xim1 = xi;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1-_a1*yim1-_a2*yim2;
      y[i] += yi;
      yim2 = yim1;
      yim1 = yi;
      xim1 = xi;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim2 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b1*xim1+_b2*xim2-_a1*yim1-_a2*yim2;
      y[i] += yi;
      yim2 = yim1;
      yim1 = yi;
      xim2 = xim1;
      xim1 = xi;
    }
  }

  // General case.
  else { 
    float yim2 = 0.0f;
    float yim1 = 0.0f;
    float xim2 = 0.0f;
    float xim1 = 0.0f;
    for (int i=0; i<n; ++i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xim1+_b2*xim2-_a1*yim1-_a2*yim2;
      y[i] += yi;
      yim2 = yim1;
      yim1 = yi;
      xim2 = xim1;
      xim1 = xi;
    }
  }
}

/**
  * Applies this filter in the reverse direction, accumulating the output. 
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must have equal
  * lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulateReverse(const vector<float>& x, vector<float>& y) {
  checkArrays(x,y);
  int n = y.size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    float yip1 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi-_a1*yip1;
      y[i] += yi;
      yip1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    float xip1 = 0.0f;
    float yip1 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1-_a1*yip1;
      y[i] += yi;
      yip1 = yi;
      xip1 = xi;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    float xip1 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1-_a1*yip1-_a2*yip2;
      y[i] += yi;
      yip2 = yip1;
      yip1 = yi;
      xip1 = xi;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    float xip1 = 0.0f;
    float xip2 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b1*xip1+_b2*xip2-_a1*yip1-_a2*yip2;
      y[i] += yi;
      yip2 = yip1;
      yip1 = yi;
      xip2 = xip1;
      xip1 = xi;
    }
  }

  // General case.
  else { 
    float xip1 = 0.0f;
    float xip2 = 0.0f;
    float yip1 = 0.0f;
    float yip2 = 0.0f;
    for (int i=n-1; i>=0; --i) {
      float xi = x[i];
      float yi = _b0*xi+_b1*xip1+_b2*xip2-_a1*yip1-_a2*yip2;
      y[i] += yi;
      yip2 = yip1;
      yip1 = yi;
      xip2 = xip1;
      xip1 = xi;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// 2-D

/**
  * Applies this filter in 1st dimension in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply1Forward(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  for (int i2=0; i2<n2; ++i2) {
    applyForward(x[i2],y[i2]);
  }
}

/**
  * Applies this filter in 1st dimension in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply1Reverse(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  for (int i2=0; i2<n2; ++i2) {
    applyReverse(x[i2],y[i2]);
  }
}

/**
  * Applies this filter in 2nd dimension in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply2Forward(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  int n1 = y[0].size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    vector<float> yim1(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> xi = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        yi[i1] = _b0*xi[i1]-
                            _a1*yim1[i1];
      }
      yim1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    vector<float> yim1(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]-
                            _a1*yim1[i1];
      }
      yim1 = yi;
      vector<float> xt = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]-
                            _a1*yim1[i1]-_a2*yim2[i1];
      }
      yim2 = yim1;
      yim1 = yi;
      vector<float> xt = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> xim2(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b1*xim1[i1]+_b2*xim2[i1]-
                  _a1*yim1[i1]-_a2*yim2[i1];
      }
      yim2 = yim1;
      yim1 = yi;
      vector<float> xt = xim2;
      xim2 = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // General case.
  else {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> xim2(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]+_b2*xim2[i1]-
                            _a1*yim1[i1]-_a2*yim2[i1];
      }
      yim2 = yim1;
      yim1 = yi;
      vector<float> xt = xim2;
      xim2 = xim1;
      xim1 = xi;
      xi = xt;
    }
  }
}

/**
  * Applies this filter in 2nd dimension in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply2Reverse(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  int n1 = y[0].size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    vector<float> yip1(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> xi = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        yi[i1] = _b0*xi[i1]-
                            _a1*yip1[i1];
      }
      yip1 = yi;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    vector<float> yip1(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]-
                            _a1*yip1[i1];
      }
      yip1 = yi;
      vector<float> xt = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]-
                            _a1*yip1[i1]-_a2*yip2[i1];
      }
      yip2 = yip1;
      yip1 = yi;
      vector<float> xt = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> xip2(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b1*xip1[i1]+_b2*xip2[i1]-
                  _a1*yip1[i1]-_a2*yip2[i1];
      }
      yip2 = yip1;
      yip1 = yi;
      vector<float> xt = xip2;
      xip2 = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // General case.
  else {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> xip2(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> yi = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]+_b2*xip2[i1]-
                            _a1*yip1[i1]-_a2*yip2[i1];
      }
      yip2 = yip1;
      yip1 = yi;
      vector<float> xt = xip2;
      xip2 = xip1;
      xip1 = xi;
      xi = xt;
    }
  }
}

/**
  * Accumulates output in 1st dimension in the forward direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate1Forward(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  for (int i2=0; i2<n2; ++i2) {
    accumulateForward(x[i2],y[i2]);
  }
}

/**
  * Accumulates output in 1st dimension in the reverse direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate1Reverse(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  for (int i2=0; i2<n2; ++i2) {
    accumulateReverse(x[i2],y[i2]);
  }
}

/**
  * Accumulates output in 2nd dimension in the forward direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate2Forward(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  int n1 = y[0].size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    vector<float> yim1(n1);
    vector<float> yi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> xi = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        yi[i1] = _b0*xi[i1]-
                            _a1*yim1[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yim1;
      yim1 = yi;
      yi = yt;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    vector<float> yim1(n1);
    vector<float> yi(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]-
                            _a1*yim1[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yim1;
      yim1 = yi;
      yi = yt;
      vector<float> xt = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> yi(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]-
                            _a1*yim1[i1]-_a2*yim2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yim2;
      yim2 = yim1;
      yim1 = yi;
      yi = yt;
      vector<float> xt = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> yi(n1);
    vector<float> xim2(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b1*xim1[i1]+_b2*xim2[i1]-
                  _a1*yim1[i1]-_a2*yim2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yim2;
      yim2 = yim1;
      yim1 = yi;
      yi = yt;
      vector<float> xt = xim2;
      xim2 = xim1;
      xim1 = xi;
      xi = xt;
    }
  }

  // General case.
  else {
    vector<float> yim2(n1);
    vector<float> yim1(n1);
    vector<float> yi(n1);
    vector<float> xim2(n1);
    vector<float> xim1(n1);
    vector<float> xi(n1);
    for (int i2=0; i2<n2; ++i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xim1[i1]+_b2*xim2[i1]-
                            _a1*yim1[i1]-_a2*yim2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yim2;
      yim2 = yim1;
      yim1 = yi;
      yi = yt;
      vector<float> xt = xim2;
      xim2 = xim1;
      xim1 = xi;
      xi = xt;
    }
  }
}

/**
/**
  * Accumulates output in 2nd dimension in the reverse direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate2Reverse(const vector< vector<float> >& x, vector< vector<float> >& y) {
  checkArrays(x,y);
  int n2 = y.size();
  int n1 = y[0].size();

  // Special case b1 = b2 = a2 = 0.
  if (_b1==0.0f && _b2==0.0f && _a2==0.0f) {
    vector<float> yip1(n1);
    vector<float> yi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> xi = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        yi[i1] = _b0*xi[i1]-
                            _a1*yip1[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yip1;
      yip1 = yi;
      yi = yt;
    }
  }

  // Special case b2 = a2 = 0.
  else if (_b2==0.0f && _a2==0.0f) {
    vector<float> yip1(n1);
    vector<float> yi(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]-
                            _a1*yip1[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yip1;
      yip1 = yi;
      yi = yt;
      vector<float> xt = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // Special case b2 = 0.
  else if (_b2==0.0f) {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> yi(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]-
                            _a1*yip1[i1]-_a2*yip2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yip2;
      yip2 = yip1;
      yip1 = yi;
      yi = yt;
      vector<float> xt = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // Special case b0 = 0.
  else if (_b0==0.0f) {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> yi(n1);
    vector<float> xip2(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b1*xip1[i1]+_b2*xip2[i1]-
                  _a1*yip1[i1]-_a2*yip2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yip2;
      yip2 = yip1;
      yip1 = yi;
      yi = yt;
      vector<float> xt = xip2;
      xip2 = xip1;
      xip1 = xi;
      xi = xt;
    }
  }

  // General case.
  else {
    vector<float> yip2(n1);
    vector<float> yip1(n1);
    vector<float> yi(n1);
    vector<float> xip2(n1);
    vector<float> xip1(n1);
    vector<float> xi(n1);
    for (int i2=n2-1; i2>=0; --i2) {
      vector<float> x2 = x[i2];
      vector<float> y2 = y[i2];
      for (int i1=0; i1<n1; ++i1) {
        xi[i1] = x2[i1];
        yi[i1] = _b0*xi[i1]+_b1*xip1[i1]+_b2*xip2[i1]-
                            _a1*yip1[i1]-_a2*yip2[i1];
        y2[i1] += yi[i1];
      }
      vector<float> yt = yip2;
      yip2 = yip1;
      yip1 = yi;
      yi = yt;
      vector<float> xt = xip2;
      xip2 = xip1;
      xip1 = xi;
      xi = xt;
    }
  }
}

///////////////////////////////////////////////////////////////////////////
// 3-D

/**
  * Applies this filter in 1st dimension in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply1Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    apply1Forward(x[i3],y[i3]);
  }
}

/**
  * Applies this filter in 1st dimension in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply1Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    apply1Reverse(x[i3],y[i3]);
  }
}

/**
  * Applies this filter in 2nd dimension in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply2Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    apply2Forward(x[i3],y[i3]);
  }
}

/**
  * Applies this filter in 2nd dimension in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply2Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    apply2Reverse(x[i3],y[i3]);
  }
}

/**
  * Applies this filter in 3rd dimension in the forward direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply3Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  int n2 = y[0].size();
  int n1 = y[0][0].size();
  vector<vector<float> > xy(n3);
  for (int i=0;i<n3; i++)
  {
    xy[i].resize(n1);
  }
  for (int i2=0; i2<n2; ++i2) {
    get2(i2,x,xy);
    apply2Forward(xy,xy);
    set2(i2,xy,y);
  }
}

/**
  * Applies this filter in 3rd dimension in the reverse direction. 
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::apply3Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  int n2 = y[0].size();
  int n1 = y[0][0].size();
  vector<vector<float> > xy(n3);
  for (int i=0;i<n3; i++)
  {
    xy[i].resize(n1);
  }
  for (int i2=0; i2<n2; ++i2) {
    get2(i2,x,xy);
    apply2Reverse(xy,xy);
    set2(i2,xy,y);
  }
}

/**
  * Accumulates output in 1st dimension in the forward direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate1Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    accumulate1Forward(x[i3],y[i3]);
  }
}

/**
  * Accumulates output in 1st dimension in the reverse direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate1Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    accumulate1Reverse(x[i3],y[i3]);
  }
}

/**
  * Accumulates output in 2nd dimension in the forward direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate2Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    accumulate2Forward(x[i3],y[i3]);
  }
}

/**
  * Accumulates output in 2nd dimension in the reverse direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate2Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  for (int i3=0; i3<n3; ++i3) {
    accumulate2Reverse(x[i3],y[i3]);
  }
}

/**
  * Accumulates output in 3rd dimension in the forward direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate3Forward(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  int n2 = y[0].size();
  int n1 = y[0][0].size();
  vector<vector<float> > xy(n3);
  for (int i=0;i<n3; i++)
  {
    xy[i].resize(n1);
  }
  for (int i2=0; i2<n2; ++i2) {
    get2(i2,x,xy);
    apply2Forward(xy,xy);
    acc2(i2,xy,y);
  }
}

/**
  * Accumulates output in 3rd dimension in the reverse direction.
  * This method filters the input, and adds the result to the output; it
  * is most useful when implementing parallel forms of recursive filters.
  * <p>
  * Input and output arrays may be the same array, but must be
  * regular and have equal lengths.
  * @param x the input array.
  * @param y the output array.
  */
void Recursive2ndOrderFilter::accumulate3Reverse(const vector<vector<vector<float> > >& x, vector<vector<vector<float> > >& y) {
  checkArrays(x,y);
  int n3 = y.size();
  int n2 = y[0].size();
  int n1 = y[0][0].size();
  vector<vector<float> > xy(n3);
  for (int i=0;i<n3; i++)
  {
    xy[i].resize(n1);
  }
  for (int i2=0; i2<n2; ++i2) {
    get2(i2,x,xy);
    apply2Reverse(xy,xy);
    acc2(i2,xy,y);
  }
}

void Recursive2ndOrderFilter::checkArrays(const vector<float>& x, const vector<float>& y) 
{
  assert(x.size()== y.size());
  assert(&x != &y);
}

void Recursive2ndOrderFilter::checkArrays(const vector< vector<float> >& x, const vector< vector<float> >& y) 
{
  assert(x.size() ==y.size());
  assert(&x != &y);
}

void Recursive2ndOrderFilter::checkArrays(const vector<vector<vector<float> > >& x, const vector<vector<vector<float> > >& y)
{
  assert(x.size() ==y.size());
  assert(&x != &y);
}

void Recursive2ndOrderFilter::get2(int i2, const vector<vector<vector<float> > >& x, vector< vector<float> >& x2) {
  int n3 = x2.size();
  int n1 = x2[0].size();
  for (int i3=0; i3<n3; ++i3) {
    vector<float> x32 = x[i3][i2];
    vector<float> x23 = x2[i3];
    for (int i1=0; i1<n1; ++i1) {
      x23[i1] = x32[i1];
    }
  }
}

void Recursive2ndOrderFilter::set2(int i2, const vector< vector<float> >& x2, vector<vector<vector<float> > >& x) {
  int n3 = x2.size();
  int n1 = x2[0].size();
  for (int i3=0; i3<n3; ++i3) {
    vector<float> x32 = x[i3][i2];
    vector<float> x23 = x2[i3];
    for (int i1=0; i1<n1; ++i1) {
      x32[i1] = x23[i1];
    }
  }
}

void Recursive2ndOrderFilter::acc2(int i2, const vector< vector<float> >& x2, vector<vector<vector<float> > >& x) {
  int n3 = x2.size();
  int n1 = x2[0].size();
  for (int i3=0; i3<n3; ++i3) {
    vector<float> x32 = x[i3][i2];
    vector<float> x23 = x2[i3];
    for (int i1=0; i1<n1; ++i1) {
      x32[i1] += x23[i1];
    }
  }
}
