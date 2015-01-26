Enter file contents here/*! \file LocalCorrelationFilter.cxx */

/*
C++ port from Java by Tiancheng Song   
Please refer to the header file for complete documentation and licensing information. 
*/
/******************************************************************************
 * Member Functions:
 *       LocalCorrelationFilter::correlate
 *       LocalCorrelationFilter::updateNormalize
 *       LocalCorrelationFilter::shift
 *       LocalCorrelationFilter::shift1
 *       LocalCorrelationFilter::shift2
 *       LocalCorrelationFilter::shift3
 *       LocalCorrelationFilter::checkDimension
 *       LocalCorrelationFilter::sub
 *       LocalCorrelationFilter::sqrtArray
 *       LocalCorrelationFilter::div
 *       LocalCorrelationFilter::multiplyArray
 *       LocalCorrelationFilter::apply2Reverse
 *       LocalCorrelationFilter::accumulate1Forward
 *       LocalCorrelationFilter::accumulate1Reverse
 *       LocalCorrelationFilter::accumulate2Forward
 *       LocalCorrelationFilter::accumulate2Reverse
 *       LocalCorrelationFilter::apply3Forward
 *       LocalCorrelationFilter::apply3Reverse
 *       GaussianFilter::apply
 *       GaussianFilter::apply1
 *       GaussianFilter::apply2
 *       GaussianFilter::apply3
 *       GaussianFilter::GaussianFilter
 ******************************************************************************
*/
#include "LocalCorrelationFilter.h"
#include <algorithm>
#include <assert.h>
#include <math.h>

#define PI 3.1415926

  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has the
   * same half-width for all dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma the correlation window half-width; must not be less than 1.
   */
  const float LocalCorrelationFilter::S1 =  0.6157280f;
  const float LocalCorrelationFilter::S2 = -0.1558022f;
  const float LocalCorrelationFilter::S3 =  0.0509014f;
  const float LocalCorrelationFilter::S4 = -0.0115417f;
  const float LocalCorrelationFilter::S[8] = {S4,S3,S2,S1,S1,S2,S3,S4};

  LocalCorrelationFilter::LocalCorrelationFilter(Type type, Window window, double sigma) 
  {
    assert(sigma>=1.0);    
    _type = type;
    _window = window;
    _sigma1 = sigma;
    _sigma2 = sigma;
    _sigma3 = sigma;
    if (window == GAUSSIAN) {
      _f1 = new GaussianFilter(sigma);
      _f2 = new GaussianFilter(sigma);
      _f3 = new GaussianFilter(sigma);
    }
#if 0 
    else {
      _f1 = new RectangleFilter(sigma);
      _f2 = new RectangleFilter(sigma);
      _f3 = new RectangleFilter(sigma);
    }
#endif
  }
  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has half-width 
   * sigma1 for the 1st dimension and half-width sigma2 for 2nd and higher 
   * dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd and higher 
   *  dimensions; must not be less than 1.
   */
  LocalCorrelationFilter::LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2) 
  {
    assert(sigma1>=1.0);
    assert(sigma2>=1.0);   
    _type = type;
    _window = window;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _sigma3 = sigma2;
    if (window== GAUSSIAN) {
      _f1 = new GaussianFilter(sigma1);
      _f2 = new GaussianFilter(sigma2);
      _f3 = new GaussianFilter(sigma2);
    } 
#if 0 
    else {
      _f1 = new RectangleFilter(sigma1);
      _f2 = new RectangleFilter(sigma2);
      _f3 = new RectangleFilter(sigma3);
    }
#endif
  }
  /**
   * Construct a correlation filter with specified parameters.
   * When applied to multi-dimensional arrays, the filter has half-width 
   * sigma1 for the 1st dimension, half-width sigma2 for the 2nd dimension,
   * and half-width sigma3 for 3rd and higher dimensions.
   * @param type the correlation type.
   * @param window the correlation window.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd dimension;
   *  must not be less than 1.
   * @param sigma3 correlation window half-width for 3rd and higher 
   * dimensions; must not be less than 1.
   */
  LocalCorrelationFilter::LocalCorrelationFilter(
    Type type, Window window, double sigma1, double sigma2, double sigma3) 
  {
    assert(sigma1>=1.0);
    assert(sigma2>=1.0);
    assert(sigma3>=1.0);
    _type = type;
    _window = window;
    _sigma1 = sigma1;
    _sigma2 = sigma2;
    _sigma3 = sigma3;
    if (window == GAUSSIAN) {
      _f1 = new GaussianFilter(sigma1);
      _f2 = new GaussianFilter(sigma2);
      _f3 = new GaussianFilter(sigma3);
    } 
#if 0 
    else {
      _f1 = new RectangleFilter(sigma1);
      _f2 = new RectangleFilter(sigma2);
      _f3 = new RectangleFilter(sigma3);
    }
#endif
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  void LocalCorrelationFilter::setInputs(const vector<float>& f, const vector<float>& g) {
    if (&f==NULL || &g==NULL) {
      _dimension = _n1 = _n2 = _n3 = 0;
    } else {
      assert(f.size()==g.size());
      _dimension = 1;
      _n1 = f.size();
      _n2 = _n3 = 0;
      _f.resize(1,vector<vector<float> >(1,vector<float>(f.size()))); 
      _g.resize(1,vector<vector<float> >(1,vector<float>(g.size()))); 
      _f[0][0] = f;
      _g[0][0] = g;
    }
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  void LocalCorrelationFilter::setInputs(const vector<vector<float> >& f, const vector<vector<float> >& g) {
    if (&f==NULL || &g==NULL) {
      _dimension = _n1 = _n2 = _n3 = 0;
    } else {
      assert(f[0].size()==g[0].size());
      assert(f.size()==g.size());
      _dimension = 2;
      _n1 = f[0].size();
      _n2 = f.size();
      _n3 = 0;
      _f.resize(1,vector<vector<float> >(f.size(),vector<float>(f[0].size()))); 
      _g.resize(1,vector<vector<float> >(g.size(),vector<float>(g[0].size()))); 
      _f[0] = f;
      _g[0] = g;
    }
  }

  /**
   * Sets the input arrays to be cross-correlated.
   * The input arrays f and g can be the same array.
   * @param f the input array f; by reference, not copied.
   * @param g the input array g; by reference, not copied.
   */
  void LocalCorrelationFilter::setInputs(const vector<vector<vector<float> > >& f, const vector<vector<vector<float> > >& g) {
    if (&f==NULL || &g==NULL) {
      _dimension = _n1 = _n2 = _n3 = 0;
    } else {
      assert(
        f[0][0].size()==g[0][0].size());
      assert(f[0].size()==g[0].size());
      assert(f.size()==g.size());
      _dimension = 3;
      _n1 = f[0][0].size();
      _n2 = f[0].size();
      _n3 = f.size();
      _f = f;
      _g = g;
    }
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag the correlation lag.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  void LocalCorrelationFilter::correlate(int lag, vector<float>& c) {
    checkDimensions(c);
    correlate(lag,_f[0][0],_g[0][0],c);
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag1 the lag in the 1st dimension.
   * @param lag2 the lag in the 2nd dimension.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  void LocalCorrelationFilter::correlate(int lag1, int lag2, vector<vector<float> >& c) {
    checkDimensions(c);
    correlate(lag1,lag2,_f[0],_g[0],c);
  }

  /**
   * Correlates the current inputs for the specified lag.
   * @param lag1 the lag in the 1st dimension.
   * @param lag2 the lag in the 2nd dimension.
   * @param lag3 the lag in the 3rd dimension.
   * @param c the output array; cannot be the same as inputs f or g.
   */
  void LocalCorrelationFilter::correlate(int lag1, int lag2, int lag3, vector<vector<vector<float> > >& c) {
    checkDimensions(c);
    correlate(lag1,lag2,lag3,_f,_g,c);
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag the lag.
   * @param c the cross-correlation to be modified.
   */
  void LocalCorrelationFilter::normalize(int lag, vector<float>& c) {
    checkDimensions(c);
    if (_s.size()==0)
      updateNormalize();
    int n1 = _n1;
    int l1 = lag;
    if (_type==SIMPLE) {
      vector<float>& sf = _s[0][0][0];
      vector<float>& sg = _s[1][0][0];
      int i1min = std::max(0,-l1);
      int i1max = std::min(n1,n1-l1);
      for (int i1=0; i1<i1min; ++i1) {
        c[i1] *= sf[i1]*sg[0];
      }
      for (int i1=i1min; i1<i1max; ++i1) {
        c[i1] *= sf[i1]*sg[i1+l1];
      }
      for (int i1=i1max; i1<n1; ++i1) {
        c[i1] *= sf[i1]*sg[n1-1];
      }
    } else if (_type==SYMMETRIC) {
      vector<float>& s = _s[0][0][0];
      for (int i1=0; i1<n1; ++i1) {
        c[i1] *= s[i1];
      }
    }
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag1 the lag.
   * @param c the cross-correlation to be modified.
   */
  void LocalCorrelationFilter::normalize(int lag1, int lag2, vector<vector<float> >& c) {
    checkDimensions(c);
    if (_s.size()==0)
      updateNormalize();
    int n1 = _n1;
    int n2 = _n2;
    int l1 = lag1;
    int l2 = lag2;
    if (_type==SIMPLE) {
      vector<vector<float> >& sf = _s[0][0];
      vector<vector<float> >& sg = _s[1][0];
      int i1min = std::max(0,-l1);
      int i1max = std::min(n1,n1-l1);
      for (int i2=0; i2<n2; ++i2) {
        vector<float>& c2 = c[i2];
        vector<float>& sf2 = sf[i2];
        vector<float>& sg2 = sg[std::max(0,std::min(n2-1,i2+l2))];
        for (int i1=0; i1<i1min; ++i1) {
          c2[i1] *= sf2[i1]*sg2[0];
        }
        for (int i1=i1min; i1<i1max; ++i1) {
          c2[i1] *= sf2[i1]*sg2[i1+l1];
        }
        for (int i1=i1max; i1<n1; ++i1) {
          c2[i1] *= sf2[i1]*sg2[n1-1];
        }
      }
    } else if (_type==SYMMETRIC) {
      vector<vector<float> >& s = _s[0][0];
      for (int i2=0; i2<n2; ++i2) {
        vector<float>& c2 = c[i2];
        vector<float>& s2 = s[i2];
        for (int i1=0; i1<n1; ++i1) {
          c2[i1] *= s2[i1];
        }
      }
    }
  }

  /**
   * Normalizes the cross-correlation for a specified lag.
   * @param lag1 the lag in the 1st dimension.
   * @param lag2 the lag in the 2nd dimension.
   * @param lag3 the lag in the 3rd dimension.
   * @param c the cross-correlation to be modified.
   */
  void LocalCorrelationFilter::normalize(int lag1, int lag2, int lag3, vector<vector<vector<float> > >& c) {
    checkDimensions(c);
    if (_s.size()==0)
      updateNormalize();
    int n1 = _n1;
    int n2 = _n2;
    int n3 = _n3;
    int l1 = lag1;
    int l2 = lag2;
    int l3 = lag3;
    if (_type==SIMPLE) {
      vector<vector<vector<float> > >& sf = _s[0];
      vector<vector<vector<float> > >& sg = _s[1];
      int i1min = std::max(0,-l1);
      int i1max = std::min(n1,n1-l1);
      for (int i3=0; i3<n3; ++i3) {
        vector<vector<float> >& c3 = c[i3];
        vector<vector<float> >& sf3 = sf[i3];
        vector<vector<float> >& sg3 = sg[std::max(0,std::min(n3-1,i3+l3))];
        for (int i2=0; i2<n2; ++i2) {
          vector<float>& c32 = c3[i2];
          vector<float>& sf32 = sf3[i2];
          vector<float>& sg32 = sg3[std::max(0,std::min(n2-1,i2+l2))];
          for (int i1=0; i1<i1min; ++i1) {
            c32[i1] *= sf32[i1]*sg32[0];
          }
          for (int i1=i1min; i1<i1max; ++i1) {
            c32[i1] *= sf32[i1]*sg32[i1+l1];
          }
          for (int i1=i1max; i1<n1; ++i1) {
            c32[i1] *= sf32[i1]*sg32[n1-1];
          }
        }
      }
    } else if (_type==SYMMETRIC) {
      vector<vector<vector<float> > >& s = _s[0];
      for (int i3=0; i3<n3; ++i3) {
        vector<vector<float> >& c3 = c[i3];
        vector<vector<float> >& s3 = s[i3];
        for (int i2=0; i2<n2; ++i2) {
          vector<float>& c32 = c3[i2];
          vector<float>& s32 = s3[i2];
          for (int i1=0; i1<n1; ++i1) {
            c32[i1] *= s32[i1];
          }
        }
      }
    }
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
  vector<float> LocalCorrelationFilter::unbias(const vector<float>& f) {
    int n1 = f.size();
    vector<float>t(n1);
    _f1->apply(f,t);
    sub(f,t,t);
    return t;
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
   vector<vector<float> > LocalCorrelationFilter::unbias(vector<vector<float> >& f) {
    int n1 = f[0].size();
    int n2 = f.size();
    vector<vector<float> > t(n2,vector<float >(n1)); 
    _f1->apply1(f,t);
    _f2->apply2(t,t);
    sub(f,t,t);
    return t;
  }

  /** 
   * Removes bias by subtracting local means from the specified array.
   * @param f the input array.
   * @return the output array, with bias subtracted.
   */
   vector<vector<vector<float> > > LocalCorrelationFilter::unbias(const vector<vector<vector<float> > >& f) {
    int n1 = f[0][0].size();
    int n2 = f[0].size();
    int n3 = f.size();
    vector<vector<vector<float> > > t(n3,vector<vector<float> >(n2,vector <float>(n1,0)));
    _f1->apply1(f,t);
    _f2->apply2(t,t);
    _f3->apply3(t,t);
    sub(f,t,t);
    return t;
  }

  void LocalCorrelationFilter::correlate(int lag, const vector<float>& f, const vector<float>& g, vector<float>& c) {
    assert(&f!=&c);
    assert(&g!=&c);
    int n1 = f.size();
    int l1 = lag;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;

    // Shifted lags for symmetric correlation.
    if (_type==SYMMETRIC) {
      // Examples of symmetric lags:
      // lag  ...  -2  -1   0   1   2  ...
      // l1f  ...  -1  -1   0   0   1  ...
      // l1g  ...  -1   0   0   1   1  ...
      l1f = (l1>=0)?(l1  )/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1  )/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    if (_window==GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==SYMMETRIC) {
      if (_window==GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
      } else {
        scale1 *= std::max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
      }
    }
    float scale = (float)scale1;
    // Correlation product.
    vector<float> h(n1);
    //int i1min = std::max(0,l1f,-l1g);
    int i1min = std::max(0, std::max(l1f,-l1g));
    //int i1max = std::min(n1,n1+l1f,n1-l1g);
    int i1max = std::min(n1,std::min(n1+l1f,n1-l1g));
    for (int i1=i1min; i1<i1max; ++i1) {
      h[i1] = scale*f[i1-l1f]*g[i1+l1g];
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==GAUSSIAN && _type==SYMMETRIC) {
      if (l1f!=l1g) {
        shift(h,c);
        //copy(c,h);
        h = c;
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    //GaussianFilter f1(*_f1);
    GaussianFilter f1(*_f1);
#if 0
    if (_window==RECTANGLE && _type==SYMMETRIC)
      f1 = new RectangleFilter(_sigma1,l1);
#endif
    f1.apply(h,c);
  }

  void LocalCorrelationFilter::correlate(
    int lag1, int lag2, const vector<vector<float> >& f, const vector<vector<float> >& g, vector<vector<float> >& c) 
  {
    assert(&f!=&c);
    assert(&g!=&c);
    int n1 = f[0].size();
    int n2 = f.size();
    int l1 = lag1;
    int l2 = lag2;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;
    int l2f = 0;
    int l2g = l2;

    // Shifted lags for symmetric correlation.
    if (_type==SYMMETRIC) {
      l1f = (l1>=0)?(l1  )/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1  )/2;
      l2f = (l2>=0)?(l2  )/2:(l2-1)/2;
      l2g = (l2>=0)?(l2+1)/2:(l2  )/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    double scale2 = 1.0;
    if (_window==GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
      scale2 *= sqrt(2.0*PI)*_sigma2;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
      scale2 *= 1.0+2.0*_sigma2;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==SYMMETRIC) {
      if (_window==GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
        scale2 *= exp((-0.125*l2*l2)/(_sigma2*_sigma2));
      } else {
        scale1 *= std::max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
        scale2 *= std::max(0.0,1.0+2.0*_sigma2-abs(l2))/(1.0+2.0*_sigma2);
      }
    }
    float scale = (float)(scale1*scale2);

    // Correlation product.
    vector<vector<float> > h(n2,vector<float >(n1)); 
    int i1min = std::max(0, std::max(l1f,-l1g));
    int i1max = std::min(n1,std::min(n1+l1f,n1-l1g));
    int i2min = std::max(0,std::max(l2f,-l2g));
    int i2max = std::min(n2,std::min(n2+l2f,n2-l2g));

    for (int i2=i2min; i2<i2max; ++i2) {
      const vector<float>& f2 = f[i2-l2f];
      const vector<float>& g2 = g[i2+l2g];
      vector<float>& h2 = h[i2];
      for (int i1=i1min; i1<i1max; ++i1) {
        h2[i1] = scale*f2[i1-l1f]*g2[i1+l1g];
      }
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==GAUSSIAN && _type==SYMMETRIC) {
      if (l1f!=l1g) {
        shift1(h,c);
        h = c;
      }
      if (l2f!=l2g) {
        shift2(h,c);        
        h = c;
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    GaussianFilter f1(*_f1);
    GaussianFilter f2(*_f2);
#if 0 
    if (_window==RECTANGLE && _type==SYMMETRIC) {
      f1 = new RectangleFilter(_sigma1,l1);
      f2 = new RectangleFilter(_sigma2,l2);
    }
#endif
    f1.apply1(h,c);
    h = c;
    f2.apply2(h,c);
  }

  void LocalCorrelationFilter::correlate(
    int lag1, int lag2, int lag3, const vector<vector<vector<float> > >& f, const vector<vector<vector<float> > >& g, vector<vector<vector<float> > >& c) 
  {
    assert(&f!=&c);
    assert(&g!=&c);
    int n1 = f[0][0].size();
    int n2 = f[0].size();
    int n3 = f.size();
    int l1 = lag1;
    int l2 = lag2;
    int l3 = lag3;

    // Conventional lags for f and g for simple correlation.
    int l1f = 0;
    int l1g = l1;
    int l2f = 0;
    int l2g = l2;
    int l3f = 0;
    int l3g = l3;

    // Shifted lags for symmetric correlation.
    if (_type==SYMMETRIC) {
      l1f = (l1>=0)?(l1  )/2:(l1-1)/2;
      l1g = (l1>=0)?(l1+1)/2:(l1  )/2;
      l2f = (l2>=0)?(l2  )/2:(l2-1)/2;
      l2g = (l2>=0)?(l2+1)/2:(l2  )/2;
      l3f = (l3>=0)?(l3  )/2:(l3-1)/2;
      l3g = (l3>=0)?(l3+1)/2:(l3  )/2;
    }

    // Scale factor so that center of window = 1.
    double scale1 = 1.0;
    double scale2 = 1.0;
    double scale3 = 1.0;
    if (_window==GAUSSIAN) {
      scale1 *= sqrt(2.0*PI)*_sigma1;
      scale2 *= sqrt(2.0*PI)*_sigma2;
      scale3 *= sqrt(2.0*PI)*_sigma3;
    } else {
      scale1 *= 1.0+2.0*_sigma1;
      scale2 *= 1.0+2.0*_sigma2;
      scale3 *= 1.0+2.0*_sigma3;
    }

    // If symmetric correlation, need extra lag-dependent scaling.
    // This scaling accounts for the separation (by lag samples) of 
    // the two windows implicitly applied to f and g. The filter we
    // apply below to the correlation product h is the product of 
    // those two windows.
    if (_type==SYMMETRIC) {
      if (_window==GAUSSIAN) {
        scale1 *= exp((-0.125*l1*l1)/(_sigma1*_sigma1));
        scale2 *= exp((-0.125*l2*l2)/(_sigma2*_sigma2));
        scale3 *= exp((-0.125*l3*l3)/(_sigma3*_sigma3));
      } else {
        scale1 *= std::max(0.0,1.0+2.0*_sigma1-abs(l1))/(1.0+2.0*_sigma1);
        scale2 *= std::max(0.0,1.0+2.0*_sigma2-abs(l2))/(1.0+2.0*_sigma2);
        scale3 *= std::max(0.0,1.0+2.0*_sigma3-abs(l3))/(1.0+2.0*_sigma3);
      }
    }
    float scale = (float)(scale1*scale2*scale3);

    // Correlation product.
    vector<vector<vector<float> > > h(n3,vector<vector<float> >(n2,vector <float>(n1))); 
    int i1min = std::max(0,std::max(l1f,-l1g));
    int i1max = std::min(n1,std::min(n1+l1f,n1-l1g));
    int i2min = std::max(0,std::max(l2f,-l2g));
    int i2max = std::min(n2,std::min(n2+l2f,n2-l2g));
    int i3min = std::max(0,std::max(l3f,-l3g));
    int i3max = std::min(n3,std::min(n3+l3f,n3-l3g));

    for (int i3=i3min; i3<i3max; ++i3) {
      const vector<vector<float> >& f3 = f[i3-l3f];
      const vector<vector<float> >& g3 = g[i3+l3g];
      vector<vector<float> >& h3 = h[i3];
      for (int i2=i2min; i2<i2max; ++i2) {
        const vector<float>& f32 = f3[i2-l2f];
        const vector<float>& g32 = g3[i2+l2g];
        vector<float>& h32 = h3[i2];
        for (int i1=i1min; i1<i1max; ++i1) {
          h32[i1] = scale*f32[i1-l1f]*g32[i1+l1g];
        }
      }
    }

    // If Gaussian and symmetric and odd lag, delay (shift) by 1/2 sample.
    if (_window==GAUSSIAN && _type==SYMMETRIC) {
      if (l1f!=l1g) {
        shift1(h,c);
        h = c;
      }
      if (l2f!=l2g) {
        shift2(h,c);   
        h = c;
      }
      if (l3f!=l3g) {
        shift3(h,c);   
        h = c;
      }
    }

    // Filter correlation product with window. For symmetric correlations
    // with a rectangle window, the width of the product rectangle depends 
    // on the lag, so we construct a new rectangle filter for each lag.
    GaussianFilter f1(*_f1);
    GaussianFilter f2(*_f2);
    GaussianFilter f3(*_f3);
#if 0
    if (_window==RECTANGLE && _type==SYMMETRIC) {
      f1 = new RectangleFilter(_sigma1,l1);
      f2 = new RectangleFilter(_sigma2,l2);
      f3 = new RectangleFilter(_sigma3,l3);
    }
#endif
    f1.apply1(h,c);
    h = c;
    f2.apply2(h,c);
    h = c;
    f3.apply3(h,c);
  }

  void LocalCorrelationFilter::updateNormalize() {
    if (_dimension==0)
      return;
    int ns = (_type==SIMPLE)?2:1;
    int n1 = std::max(1,_n1);
    int n2 = std::max(1,_n2);
    int n3 = std::max(1,_n3);
    _s = vector<vector<vector<vector<float> > > >(ns, vector<vector<vector<float> > >(n3, vector<vector<float> >(n2,vector <float>(n1))));
    if (_type==SIMPLE) {
      if (_dimension==1) {
        vector<float>& f = _f[0][0];
        vector<float>& g = _g[0][0];
        vector<float>& sf = _s[0][0][0];
        vector<float>& sg = _s[1][0][0];
        correlate(0,f,f,sf);
        correlate(0,g,g,sg);
        sqrtArray(sf,sf);
        sqrtArray(sg,sg);
        div(1.0f,sf,sf);
        div(1.0f,sg,sg);
      } else if (_dimension==2) {
        vector<vector<float> >& f = _f[0];
        vector<vector<float> >& g = _g[0];
        vector<vector<float> >& sf = _s[0][0];
        vector<vector<float> >& sg = _s[1][0];
        correlate(0,0,f,f,sf);
        correlate(0,0,g,g,sg);
        sqrtArray(sf,sf);
        sqrtArray(sg,sg);
        div(1.0f,sf,sf);
        div(1.0f,sg,sg);
      } else {
        vector<vector<vector<float> > >&  f = _f;
        vector<vector<vector<float> > >&  g = _g;
        vector<vector<vector<float> > >&  sf = _s[0];
        vector<vector<vector<float> > >&  sg = _s[1];
        correlate(0,0,0,f,f,sf);
        correlate(0,0,0,g,g,sg);
        sqrtArray(sf,sf);
        sqrtArray(sg,sg);
        div(1.0f,sf,sf);
        div(1.0f,sg,sg);
      }
    } else {
      if (_dimension==1) {
        vector<float>& f = _f[0][0];
        vector<float>& g = _g[0][0];
        vector<float>& s = _s[0][0][0];
        vector<float>& sf = s;
        vector<float> sg(_n1); 
        correlate(0,f,f,sf);
        correlate(0,g,g,sg);
        multiplyArray(sf,sg,s);
        sqrtArray(s,s);
        div(1.0f,s,s);
      } else if (_dimension==2) {
        vector<vector<float> >& f = _f[0];
        vector<vector<float> >& g = _g[0];
        vector<vector<float> >& s = _s[0][0];
        vector<vector<float> >& sf = s;
        vector<vector<float> > sg(_n2, vector<float>(_n1));
        correlate(0,0,f,f,sf);
        correlate(0,0,g,g,sg);
        multiplyArray(sf,sg,s);
        sqrtArray(s,s);
        div(1.0f,s,s);
      } else {
        vector<vector<vector<float> > >&  f = _f;
        vector<vector<vector<float> > >&  g = _g;
        vector<vector<vector<float> > >&  s = _s[0];
        vector<vector<vector<float> > >&  sf = s;
        vector<vector<vector<float> > > sg(n3,vector<vector<float> >(n2,vector <float>(n1,0)));
        correlate(0,0,0,f,f,sf);
        correlate(0,0,0,g,g,sg);
        multiplyArray(sf,sg,s);
        sqrtArray(s,s);
        div(1.0f,s,s);
      }
    }
  }

  void LocalCorrelationFilter::shift(const vector<float>& f, vector<float>& g) {
    int n1 = f.size();
    int i1b,i1e;

    // Rolling on.
    i1b = 0;
    i1e = std::min(4,n1);
    for (int i1=i1b; i1<i1e; ++i1) {
      int ib = std::max(0,4-i1);
      int ie = std::min(8,4-i1+n1);
      g[i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        g[i1] += S[i]*f[i1+i-4];
    }

    // Middle.
    i1b = 4;
    i1e = n1-3;
    for (int i1=i1b; i1<i1e; ++i1) {
      g[i1] = S4*(f[i1-4]+f[i1+3]) +
              S3*(f[i1-3]+f[i1+2]) +
              S2*(f[i1-2]+f[i1+1]) +
              S1*(f[i1-1]+f[i1  ]);
    }

    // Rolling off.
    i1b = std::max(0,n1-3);
    i1e = n1;
    for (int i1=i1b; i1<i1e; ++i1) {
      int ib = std::max(0,4-i1);
      int ie = std::min(8,4-i1+n1);
      g[i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        g[i1] += S[i]*f[i1+i-4];
    }
  }

  void LocalCorrelationFilter::shift1(const vector<vector<float> >& f, vector<vector<float> >& g) {
    int n2 = f.size();
    for (int i2=0; i2<n2; ++i2)
      shift(f[i2],g[i2]);
  }

  void LocalCorrelationFilter::shift2(const vector<vector<float> >& f, vector<vector<float> >& g) {
    int n2 = f.size();
    int n1 = f[0].size();
    int i2b,i2e;

    // Rolling on.
    i2b = 0;
    i2e = std::min(4,n2);
    for (int i2=i2b; i2<i2e; ++i2) {
      int ib = std::max(0,4-i2);
      int ie = std::min(8,4-i2+n2);
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] += S[i]*f[i2+i-4][i1];
    }

    // Middle.
    i2b = 4;
    i2e = n2-3;
    for (int i2=i2b; i2<i2e; ++i2) {
      vector<float>& g2 = g[i2];
      const vector<float>& fm4 = f[i2-4];
      const vector<float>& fm3 = f[i2-3];
      const vector<float>& fm2 = f[i2-2];
      const vector<float>& fm1 = f[i2-1];
      const vector<float>& fp0 = f[i2  ];
      const vector<float>& fp1 = f[i2+1];
      const vector<float>& fp2 = f[i2+2];
      const vector<float>& fp3 = f[i2+3];
      for (int i1=0; i1<n1; ++i1)
        g2[i1] = S4*(fm4[i1]+fp3[i1]) +
                 S3*(fm3[i1]+fp2[i1]) +
                 S2*(fm2[i1]+fp1[i1]) +
                 S1*(fm1[i1]+fp0[i1]);
    }

    // Rolling off.
    i2b = std::max(0,n2-3);
    i2e = n2;
    for (int i2=i2b; i2<i2e; ++i2) {
      int ib = std::max(0,4-i2);
      int ie = std::min(8,4-i2+n2);
      for (int i1=0; i1<n1; ++i1)
        g[i2][i1] = 0.0f;
      for (int i=ib; i<ie; ++i)
        for (int i1=0; i1<n1; ++i1)
          g[i2][i1] += S[i]*f[i2+i-4][i1];
    }
  }

  void LocalCorrelationFilter::shift1(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g) {
    int n3 = f.size();
    for (int i3=0; i3<n3; ++i3)
      shift1(f[i3],g[i3]);
  }

  void LocalCorrelationFilter::shift2(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g) {
    int n3 = f.size();
    for (int i3=0; i3<n3; ++i3)
      shift2(f[i3],g[i3]);
  }

  void LocalCorrelationFilter::shift3(const vector<vector<vector<float> > >& f, vector<vector<vector<float> > >& g) {
    int n3 = f.size();
    int n2 = f[0].size();
    vector<vector<float> > f2(n3,vector<float >(f[0][0].size())); 
    vector<vector<float> > g2(n3,vector<float >(g[0][0].size())); 
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {
        f2[i3] = f[i3][i2];
        g2[i3] = g[i3][i2];
      }
      shift2(f2,g2);
    }
    // Set back the result
    for (int i2=0; i2<n2; ++i2) {
      for (int i3=0; i3<n3; ++i3) {    
        g[i3][i2] = g2[i3];
      }     
    }
  }

  void LocalCorrelationFilter::checkDimension(int dimension) {
    assert(_dimension==dimension);
  }

  void LocalCorrelationFilter::checkDimensions(vector<float>& c) {
    assert(_n1==c.size());
    checkDimension(1);
  }

  void LocalCorrelationFilter::checkDimensions(vector<vector<float> >& c) {
    assert(_n1==c[0].size());
    assert(_n2==c.size());
    checkDimension(2);
  }

  void LocalCorrelationFilter::checkDimensions(vector<vector<vector<float> > >& c) {
    assert(_n1==c[0][0].size());
    assert(_n2==c[0].size());
    assert(_n3==c.size());
    checkDimension(3);
  }


void LocalCorrelationFilter::sub(const vector<float>&a,const vector<float>&b, vector<float>&c)
{
  for (unsigned i= 0; i<a.size(); i++)
    c[i] = a[i]-b[i];
}
void LocalCorrelationFilter::sub(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
      c[i][j] = a[i][j]-b[i][j];
}
void LocalCorrelationFilter::sub(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
          for (unsigned k= 0; k<a[0][0].size(); k++)
            c[i][j][k] = a[i][j][k]-b[i][j][k];
}
void LocalCorrelationFilter::sqrtArray(const vector<float>& a, vector<float>& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    c[i] = sqrt(a[i]);
}
void LocalCorrelationFilter::sqrtArray(const vector<vector<float> >& a, vector<vector<float> >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
      c[i][j] = sqrt(a[i][j]);
}
void LocalCorrelationFilter::sqrtArray(const vector<vector<vector<float> > >& a, vector<vector<vector<float> > >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
          for (unsigned k= 0; k<a[0][0].size(); k++)
             c[i][j][k] = sqrt(a[i][j][k]);
}
void LocalCorrelationFilter::div(const vector<float>&a,const vector<float>&b, vector<float>&c)
{
  for (unsigned i= 0; i<a.size(); i++)
    c[i] = a[i]/b[i];
}
void LocalCorrelationFilter::div(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
      c[i][j] = a[i][j]/b[i][j];
}
void LocalCorrelationFilter::div(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
          for (unsigned k= 0; k<a[0][0].size(); k++)
            c[i][j][k] = a[i][j][k]/b[i][j][k];
}
  
void LocalCorrelationFilter::div(const float a,const vector<float>&b, vector<float>&c)
{
  for (unsigned i= 0; i<b.size(); i++)
    c[i] = a/b[i];
}
void LocalCorrelationFilter::div(const float a,const vector<vector<float> >& b, vector<vector<float> >& c)
{
  for (unsigned i= 0; i<b.size(); i++)
    for (unsigned j= 0; j<b[0].size(); j++)
      c[i][j] = a/b[i][j];
}
void LocalCorrelationFilter::div(const float a,const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c)
{
  for (unsigned i= 0; i<b.size(); i++)
    for (unsigned j= 0; j<b[0].size(); j++)
          for (unsigned k= 0; k<b[0][0].size(); k++)
            c[i][j][k] = a/b[i][j][k];
}

void LocalCorrelationFilter::multiplyArray(const vector<float>&a,const vector<float>&b, vector<float>&c)
{
  for (unsigned i= 0; i<a.size(); i++)
    c[i] = a[i]*b[i];
}

void LocalCorrelationFilter::multiplyArray(const vector<vector<float> >& a,const vector<vector<float> >& b, vector<vector<float> >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
      c[i][j] = a[i][j]*b[i][j];
}

void LocalCorrelationFilter::multiplyArray(const vector<vector<vector<float> > >& a, const vector<vector<vector<float> > >& b, vector<vector<vector<float> > >& c)
{
  for (unsigned i= 0; i<a.size(); i++)
    for (unsigned j= 0; j<a[0].size(); j++)
      for (unsigned k= 0; k<a[0][0].size(); k++)
            c[i][j][k] = a[i][j][k]*b[i][j][k];
}
