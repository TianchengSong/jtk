/*! \file Sampling.cxx */

/*
C++ port from Java by Tiancheng Song   
Please refer to the header file for complete documentation and licensing information. 
*/

/******************************************************************************
 * Member Functions:
 *       Sampling::Sampling
 *       Sampling::getCount
 *       Sampling::getDelta
 *       Sampling::getFirst
 *       Sampling::getLast
 *       Sampling::getValue
 *       Sampling::getValues
 *       Sampling::isUniform
 *       Sampling::isEquivalentTo
 *       Sampling::isCompatible
 *       Sampling::indexOf
 *       Sampling::indexOfNearest
 *       Sampling::valueOfNearest
 *       Sampling::isInBounds
 *       Sampling::isInBounds
 *       Sampling::isInBoundsExtended
 *       Sampling::getValueExtended
 *       Sampling::indexOfNearestExtended
 *       Sampling::valueOfNearestExtended
 *       Sampling::indexOfFloorExtended
 *       Sampling::normalizedDifferenceExtended
 *       Sampling::overlapWith
 *       Sampling::mergeWith
 *       Sampling::shift
 *       Sampling::prepend
 *       Sampling::append
 *       Sampling::decimate
 *       Sampling::interpolate
 ******************************************************************************
*/

#include <math.h>
#include <string>
#include <assert.h>
#include <algorithm>  
#include "Sampling.h"
using namespace std;
/**
 * Sampling of one variable.
 * <p>
 * Samplings are often used to represent independent variables for sampled 
 * functions. They describe the values at which a function is sampled. For 
 * efficiency, and to guarantee a unique mapping from sample value to 
 * function value, we restrict samplings to be strictly increasing. In other
 * words, no two samples have equal value, and sample values increase with 
 * increasing sample index.
 * <p>
 * Samplings are either uniform or non-uniform. Uniform samplings are
 * represented by a sample count n, a sampling interval d, and a first
 * sample value f. Non-uniform samplings are represented by an array of 
 * sample values.
 * <p>
 * All sample values are computed and stored in <em>double precision</em>. 
 * This double precision can be especially important in uniform samplings, 
 * where the sampling interval d and first sample value f may be used to
 * compute values for thousands of samples, in loops like this one:
 * <pre><code>
 *   int n = sampling.getCount();
 *   double d = sampling.getDelta();
 *   double f = sampling.getFirst();
 *   double v = f;
 *   for (int i=0; i&lt;n; ++i,v+=d) {
 *     // some calculation that uses the sample value v
 *   }
 * </code></pre>
 * In each iteration of the loop above, the sample value v is computed by 
 * accumulating the sampling interval d. This computation is fast, but it 
 * also yields rounding error that can grow quadratically with the number 
 * of samples n. If v were computed in single (float) precision, then this 
 * rounding error could exceed the sampling interval d for as few as 
 * n=10,000 samples.
 * <p>
 * If accumulating in double precision is insufficient, a more accurate 
 * and more costly way to compute sample values is as follows:
 * <pre><code>
 *   // ...
 *   double v = f;
 *   for (int i=0; i&lt;n; ++i,v=f+i*d) {
 *     // some calculation that uses the sample value v
 *   }
 * </code></pre>
 * With this computation of sample values, rounding errors can grow only 
 * linearly with the number of samples n.
 * <p>
 * Two samplings are considered equivalent if their sample values differ
 * by no more than the <em>sampling tolerance</em>. This tolerance may be
 * specified, as a fraction of the sampling interval, when a sampling is
 * constructed. Alternatively, a default tolerance may be used. When
 * comparing two samplings, the smaller of their tolerances is used.
 * <p>
 * A sampling is immutable. New samplings can be constructed by applying
 * various transformations (e.g., shifting) to an existing sampling, but
 * an existing sampling cannot be changed. Therefore, multiple sampled
 * functions can safely share the same sampling.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2005.03.11
 */


  /**
   * A default fraction used to test for equivalent sample values. 
   * By default, if the difference between two sample values does not 
   * exceed this fraction of the sampling interval, then those values 
   * are considered equivalent. This default is used when a tolerance 
   * is not specified explicitly when a sampling is constructed.
   */
  double Sampling::DEFAULT_TOLERANCE = 1.0e-6;

  /**
   * Constructs a uniform sampling with specified count. The sampling 
   * interval (delta) is 1.0, and the first sample value is 0.0.
   * @param n the number (count) of samples; must be positive.
   */

  void Sampling::init(int n, double d, double f, double t)
  {
    assert(n>0);
    assert(d>0.0);
    _n = n;
    _d = d;
    _f = f;
    _v.resize(0);
    _t = t;
    _td = _t*_d;
  }

  Sampling::Sampling(int n) {
    init(n, 1.0, 0.0, DEFAULT_TOLERANCE);    
  }

  /**
   * Constructs a uniform sampling with specified parameters.
   * @param n the number (count) of samples; must be positive.
   * @param d the sampling interval (delta); must be positive.
   * @param f the first sample value.
   */
  Sampling::Sampling(int n, double d, double f) {
    init(n, d, f, DEFAULT_TOLERANCE);    
  }

  /**
   * Constructs a sampling with specified parameters.
   * @param n the number (count) of samples; must be positive.
   * @param d the sampling interval (delta); must be positive.
   * @param f the first sample value.
   * @param t the sampling tolerance, expressed as fraction of delta.
   */
  Sampling::Sampling(int n, double d, double f, double t) {
    init(n, d, f, t);
  }

  /**
   * Constructs a sampling from the specified array of values. The values 
   * must be strictly increasing.
   * <p>
   * The constructed sampling may or may not be uniform, depending on the 
   * specified values and default sampling tolerance. If uniform (to within 
   * the default tolerance), then the array of values is discarded, and the 
   * sampling is represented by the count, sampling interval, and first 
   * sample value.
   * @param v the array of sampling values; must have non-zero length.
   */
  Sampling::Sampling(vector<double>& v) {
    //this(v,DEFAULT_TOLERANCE);
    assert(v.size()>0);
   // assert(isIncreasing(v));
    bool testing = true;
    if (v.size() > 1)
    {
      for (unsigned i = 0; i<v.size(); i++)
      {
        if (v[i]<v[i-1])
        {
          testing = false;
          break;
        }
      }
    }
    assert(testing);

    _n = v.size();
    _d = (_n<2)?1.0:(v[_n-1]-v[0])/(_n-1);
    _f = v[0];
    _t = DEFAULT_TOLERANCE;
    _td = _t*_d;
    bool uniform = true;
    for (int i=0; i<_n && uniform; ++i) {
      double vi = _f+i*_d;
      if (!almostEqual(v[i],vi,_td))
        uniform = false;
    }
    if (uniform) 
      _v.resize(0);
    else
      _v = v;
  }

  /**
   * Constructs a sampling from the specified array of values and tolerance.
   * The values must be strictly increasing. 
   * <p>
   * The constructed sampling may or may not be uniform, depending on the 
   * specified values and tolerance. If uniform (to within the specified
   * tolerance), then the array of values is discarded, and the sampling is 
   * represented by the count, sampling interval, and first sample value.
   * @param v the array of sampling values; must have non-zero length.
   * @param t the sampling tolerance, expressed as fraction of delta.
   */
  Sampling::Sampling(vector<double>& v, double t) {
    assert(v.size()>0);
    bool testing = true;
    if (v.size() > 1)
    {
      for (unsigned i = 0; i<v.size(); i++)
      {
        if (v[i]<v[i-1])
        {
          testing = false;
          break;
        }
      }
    }
    assert(testing);

    _n = v.size();
    _d = (_n<2)?1.0:(v[_n-1]-v[0])/(_n-1);
    _f = v[0];
    _t = t;
    _td = _t*_d;
    bool uniform = true;
    for (int i=0; i<_n && uniform; ++i) {
      double vi = _f+i*_d;
      if (!almostEqual(v[i],vi,_td))
        uniform = false;
    }
    if (uniform) 
      _v.resize(0);
    else
      _v = v;
  }

  /**
   * Gets the number of samples.
   * @return the number of samples.
   */
  int Sampling::getCount() {
    return _n;
  }

  /**
   * Gets the sampling interval. If not uniformly sampled, the sampling
   * interval is the average difference between sample values.
   * @return the sampling interval; 1.0, if fewer than two samples.
   */
  double Sampling::getDelta() {
    return _d;
  }

  /**
   * Gets the first sample value.
   * @return the first sample value.
   */
  double Sampling::getFirst() {
    return _f;
  }

  /**
   * Gets the last sample value.
   * @return the last sample value.
   */
  double Sampling::getLast() {
    //return (_v!=NULL)?_v[_n-1]:_f+(_n-1)*_d;
    return (_v.size() != 0)?_v[_n-1]:_f+(_n-1)*_d;
  }

  /**
   * Gets the sample value with specified index.
   * @param i the index.
   * @return the sample value.
   */
  double Sampling::getValue(int i) {
    assert(i<_n);
    return value(i);
  }

  /**
   * Gets the sample values. If this sampling was constructed with an array 
   * of sample values, then the returned values are equivalent (equal to 
   * within the sampling tolerance) to the values in that array.
   * @return the sample values; returned by copy, not by reference.
   */
  vector<double> Sampling::getValues() {
    vector<double> v = _v;
    if (_v.size() != 0) {
      v  = _v;
    } 
    if (v.size() == 0)
    {      
      for (int i=0; i<_n; ++i)
        v.push_back(_f+i*_d);
    }
    return v;
  }

  /**
   * Determines whether this sampling is uniform. A sampling is uniform 
   * if its values can be computed, to within the sampling tolerance, by
   * the expression v = f+i*d, for sampling indices i = 0, 1, ..., n-1.
   * Samplings with only one sample are considered to be uniform.
   * <p>
   * Note that, by this definition, samplings constructed with an array 
   * of sample values may or may not be uniform.
   * @return true, if uniform; false, otherwise.
   */
  bool Sampling::isUniform() {
    return _v.size() == 0;
  }

  /**
   * Determines whether this sampling is equivalent to the specified sampling.
   * Two samplings are equivalent if each of their sample values differs by 
   * no more than the sampling tolerance.
   * @param s the sampling to compare to this sampling.
   * @return true, if equivalent; false, otherwise.
   */
  bool Sampling::isEquivalentTo(Sampling s) {
    Sampling t(*this);
    if (t.isUniform()!=s.isUniform())
      return false;
    if (t.isUniform()) {
      if (t.getCount()!=s.getCount())
        return false;
      double tiny = tinyWith(s);
      double tf = t.getFirst();
      double tl = t.getLast();
      double sf = s.getFirst();
      double sl = s.getLast();
      return almostEqual(tf,sf,tiny) && almostEqual(tl,sl,tiny);
    } else {
      double tiny = tinyWith(s);
      for (int i=0; i<_n; ++i) {
        if (!almostEqual(_v[i],s.value(i),tiny))
          return false;
      }
      return true;
    }
  }

  /**
   * Determines whether this sampling is compatible with the specified sampling.
   * Two samplings are incompatible if their ranges of sample values overlap, 
   * but not all values in the overlapping parts are equivalent. Otherwise,
   * they are compatible.
   * @param s the sampling to compare to this sampling.
   * @return true, if compatible; false, otherwise.
   */
  bool Sampling::isCompatible(Sampling s) {
    Sampling t(*this);
    int nt = t.getCount();
    int ns = s.getCount();
    double tf = t.getFirst();
    double sf = s.getFirst();
    double tl = t.getLast();
    double sl = s.getLast();
    int it = 0;
    int is = 0;
    int jt = nt-1;
    int js = ns-1;
    if (tl<sf) {
      return true;
    } else if (sl<tf) {
      return true;
    } else {
      if (tf<sf) {
        it = t.indexOf(sf);
      } else {
        is = s.indexOf(tf);
      }
      if (it<0 || is<0)
        return false;
      if (tl<sl) {
        js = s.indexOf(tl);
      } else {
        jt = t.indexOf(sl);
      }
      if (jt<0 || js<0)
        return false;
      int mt = 1+jt-it;
      int ms = 1+js-is;
      if (mt!=ms)
        return false;
      if (!t.isUniform() || !s.isUniform()) {
        double tiny = tinyWith(s);
        for (jt=it,js=is; jt!=mt; ++jt,++js) {
          if (!almostEqual(t.value(jt),s.value(js),tiny))
            return false;
        }
      }
      return true;
    }
  }

  /**
   * Returns the index of the sample with specified value. If this 
   * sampling has a sample value that equals (to within the sampling 
   * tolerance) the specified value, then this method returns the 
   * index of that sample. Otherwise, this method returns -1.
   * @param x the value.
   * @return the index of the matching sample; -1, if none.
   */
  int Sampling::indexOf(double x) {
    int i = -1;
    if (isUniform()) {
      int j = (int)round((x-_f)/_d);
      if (0<=j && j<_n && almostEqual(x,_f+j*_d,_td))
        i = j;
    } else {
      int j = std::binary_search(_v.begin(), _v.end(), x);
      if (0<=j) {
        i = j;
      } else {
        j = -(j+1);
        if (j>0 && almostEqual(x,_v[j-1],_td)) {
          i = j-1;
        } else if (j<_n && almostEqual(x,_v[j],_td)) {
          i = j;
        }
      }
    }
    return i;
  }

  /**
   * Returns the index of the sample nearest to the specified value.
   * @param x the value.
   * @return the index of the nearest sample.
   */
  int Sampling::indexOfNearest(double x) {
    int i;
    if (isUniform()) {
      i = (int)round((x-_f)/_d);
      if (i<0)
        i = 0;
      if (i>=_n)
        i = _n-1;
    } else {
      i = std::binary_search(_v.begin(), _v.end(), x);
      if (i<0) {
        i = -(i+1);
        if (i==_n) {
          i = _n-1;
        } else if (i>0 && abs(x-_v[i-1])<abs(x-_v[i])) {
          --i;
        }
      }
    }
    return i;
  }

  /**
   * Returns the value of the sample nearest to the specified value.
   * @param x the value.
   * @return the value of the nearest sample.
   */
  double Sampling::valueOfNearest(double x) {
    return getValue(indexOfNearest(x));
  }

  /**
   * Determines whether the specified index is in the bounds of this sampling.
   * An index is in bounds if in the range [0,count-1] of the first and last
   * sample indices.
   * @param i the index.
   * @return true, if in bounds; false, otherwise.
   */
  bool Sampling::isInBounds(int i) {
    return 0<=i && i<_n; 
  }

  /**
   * Determines whether the specified value is in the bounds of this sampling.
   * A value is in bounds if in the range [first,last] defined by the first 
   * and last sample values.
   * @param x the value.
   * @return true, if in bounds; false, otherwise.
   */
  bool Sampling::isInBounds(double x) {
    return getFirst()<=x && x<=getLast();
  }

  /**
   * Determines whether the specified value is in the bounds of this sampling,
   * which is assumed to be uniform. A value is in bounds if in the range
   * [first-0.5*delta,last+0.5*delta] defined by the first and last sample 
   * values and the sampling interval delta. In effect, this method extends
   * the bounds of this sampling by one-half sample when testing the value.
   * @param x the value.
   * @return true, if in bounds; false, otherwise.
   */
  bool Sampling::isInBoundsExtended(double x) {
    assert(isUniform());
    double dhalf = 0.5*_d;
    return getFirst()-dhalf<=x && x<=getLast()+dhalf;
  }

  /**
   * Gets the value for the specified index, assuming uniform sampling.
   * The index and the returned value need not be in the bounds of this 
   * sampling, which must be uniform. That is, the specified index may be 
   * less than zero or greater than or equal to the number of samples.
   * @param i the index.
   * @return the value.
   */
  double Sampling::getValueExtended(int i) {
    assert(isUniform());
    return _f+i*_d;
  }

  /**
   * Returns the index of the sample nearest the specified value, assuming 
   * uniform sampling. The value and the returned index need not be in the 
   * bounds of this sampling, which must be uniform. Specifically, the
   * returned index may be less than zero or greater than or equal to the 
   * number of samples.
   * @param x the value.
   */
  int Sampling::indexOfNearestExtended(double x) {
    assert(isUniform());
    return (int)round((x-_f)/_d);
  }

  /**
   * Returns the value of the sample nearest to the specified value, 
   * assuming uniform sampling. The specified and returned values need 
   * not be in the bounds of this sampling, which must be uniform.
   * @param x the value.
   * @return the value of the nearest sample.
   */
  double Sampling::valueOfNearestExtended(double x) {
    return getValueExtended(indexOfNearestExtended(x));
  }

  /**
   * Returns the index of the floor of the specified value.
   * The floor is the largest sampled value less than or equal to the 
   * specified value. The value and the returned index need not be in 
   * the bounds of this sampling, which must be uniform.
   * Specifically, the returned index may be less than zero or greater 
   * than or equal to the number of samples.
   * @param x the value.
   */
  int Sampling::indexOfFloorExtended(double x) {
    double xn = (x-_f)/_d;
    int ix = (int)xn;
    return (ix<=xn)?ix:ix-1;
  }

  /**
   * Returns the normalized difference between a specified value and 
   * the sampled value for a specified index. Normalized difference is
   * the difference (specified value minus sampled value) divided by 
   * the sampling interval. The specified value and index not be in 
   * the bounds of this sampling, which must be uniform.
   * @param x the value.
   * @param i the index.
   */
  double Sampling::normalizedDifferenceExtended(double x, int i) {
    assert(isUniform());
    return (x-(_f+i*_d))/_d;
  }

  /**
   * Determines the overlap between this sampling and the specified sampling.
   * Both the specified sampling and this sampling represent a first-to-last
   * range of sample values. The overlap is a contiguous set of samples that 
   * have values that are equal, to within the minimum sampling tolerance of 
   * the two samplings. This set is represented by an array of three ints, 
   * {n,it,is}, where n is the number of overlapping samples, and it and is 
   * denote the indices of the first samples in the overlapping parts of this 
   * sampling (t) and the specified sampling (s). There exist three cases.
   * <ul><li>
   * The ranges of sample values overlap, and all values in the overlapping 
   * parts are equivalent. In this case, the two samplings are compatible; 
   * and this method returns the array {n,it,is}, as described 
   * above.
   * </li><li>
   * The ranges of sample values in the two samplings do not overlap.
   * In this case, the two samplings are compatible; and this method 
   * returns either {0,nt,0} or {0,0,ns}, depending on whether all
   * sample values in this sampling are less than or greater than
   * those in the specified sampling, respectively.
   * </li><li>
   * The ranges of sample values overlap, but not all values in the 
   * overlapping parts are equivalent. In this case, the two samplings 
   * are incompatible; and this method returns NULL.
   * </li></ul>
   * @param s the sampling to compare with this sampling.
   * @return the array {n,it,is} that describes the overlap; NULL, if 
   *  the specified sampling is incompatible with this sampling.
   */
    vector<int> Sampling::overlapWith(Sampling s) {
    Sampling t(*this);
    int nt = t.getCount();
    int ns = s.getCount();
    double tf = t.getFirst();
    double sf = s.getFirst();
    double tl = t.getLast();
    double sl = s.getLast();
    int it = 0;
    int is = 0;
    int jt = nt-1;
    int js = ns-1;
    if (tl<sf) {
      vector<int>result(3); //new int[]{0,nt,0};
      result[0] = 0;
      result[1] = nt;
      result[2] = 0;
      return result;
    } else if (sl<tf) {
      vector<int>result(3);
      result[0] = 0;
      result[1] = 0;
      result[2] = ns;
      return result;
    } else {
      if (tf<sf) {
        it = t.indexOf(sf);
      } else {
        is = s.indexOf(tf);
      }
      if (it<0 || is<0)
      {
        vector<int>result(0);
        return result;
      }

      if (tl<sl) {
        js = s.indexOf(tl);
      } else {
        jt = t.indexOf(sl);
      }
      if (jt<0 || js<0)
      {
        vector<int>result(0);
        return result;
      }
      int mt = 1+jt-it;
      int ms = 1+js-is;
      if (mt!=ms)
      {
        vector<int>result(0);
        return result;
      }       
      if (!t.isUniform() || !s.isUniform()) {
        double tiny = tinyWith(s);
        for (jt=it,js=is; jt!=mt; ++jt,++js) {
          if (!almostEqual(t.value(jt),s.value(js),tiny))
          {
            vector<int>result(0);
            return result;
          }
        }
      }
      vector<int>result(3);
      result[0] = mt;
      result[1] = it;
      result[2] = is;
      return result;
    }
  }

  /**
   * Returns the union of this sampling with the specified sampling. This
   * union is possible if and only if the two samplings are compatible.
   * <p>
   * If the two samplings do not overlap, this method does not create
   * samples within any gap that may exist between them. In other words,
   * the number of samples in the sampling returned is exactly nt+ns-n, 
   * where nt is the number of samples in this sampling, ns is the number 
   * of samples in the specified sampling, and n is the number of samples 
   * with equivalent values in any overlapping parts of the two samplings.
   * If the samplings do not overlap, then n = 0. One consequence of this
   * behavior is that the union of two uniform samplings with the same 
   * sampling interval may be non-uniform.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @see #overlapWith(Sampling)
   * @param s the sampling to merge with this sampling.
   * @return the merged sampling; NULL, if no merge is possible.
   */
  Sampling Sampling:: mergeWith(Sampling s) {
    Sampling t(*this);
    vector<int> overlap = t.overlapWith(s);
    if (overlap.size() == 0)
      return NULL;
    int n = overlap[0];
    int it = overlap[1];
    int is = overlap[2];
    int nt = t.getCount();
    int ns = s.getCount();
    int nm = nt+ns-n;
    if (n>0 && t.isUniform() && s.isUniform()) {
      double dm = t.getDelta();
      double fm = (it==0)?s.getFirst():t.getFirst();
      return  Sampling(nm,dm,fm);
    } else {
      vector<double> vm(nm); //= new double[nm];
      int jm = 0;
      for (int jt=0; jt<it; ++jt)
        vm[jm++] = t.value(jt);
      for (int js=0; js<is; ++js)
        vm[jm++] = s.value(js);
      for (int jt=it; jt<it+n; ++jt)
        vm[jm++] = t.value(jt);
      for (int jt=it+n; jt<nt; ++jt)
        vm[jm++] = t.value(jt);
      for (int js=is+n; js<ns; ++js)
        vm[jm++] = s.value(js);
      return  Sampling(vm);
    }
  }

  /**
   * Shifts this sampling.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @param s the value (shift) to add to this sampling's values.
   * @return the new sampling.
   */
  Sampling Sampling:: shift(double s) {
    if (_v.size() == 0) {
      return Sampling(_n,_d,_f+s,_t);
    } else {
      vector<double> v(_n);
      for (int i=0; i<_n; ++i)
        v[i] = _v[i]+s;
      return  Sampling(v,_t);
    }
  }

  /**
   * Prepends samples to this sampling.
   * If this sampling is not uniform, prepended sample values are computed 
   * to preserve the average difference between adjacent sample values.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @param m the number of new samples prepended to this sampling.
   * @return the new sampling.
   */
  Sampling Sampling:: prepend(int m) {
    int n = _n+m;
    double f = _f-m*_d;
    if (_v.size() == 0) {
      return Sampling(n,_d,f,_t);
    } else {
      vector<double> v (n);
      for (int i=0; i<m; ++i)
        v[i] = f+i*_d;
      for (int i=m; i<n; ++i)
        v[i] = _v[i-m];
      return  Sampling(v,_t);
    }
  }

  /**
   * Appends samples to this sampling.
   * If this sampling is not uniform, appended sample values are computed 
   * to preserve the average difference between adjacent sample values.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @param m the number of new samples appended to this sampling.
   * @return the new sampling.
   */
  Sampling Sampling:: append(int m) {
    int n = _n+m;
    if (_v.size() == 0) {
      return  Sampling(n,_d,_f,_t);
    } else {
      vector<double> v(n);
      for (int i=0; i<_n; ++i)
        v[i] = _v[i];
      for (int i=_n; i<n; ++i)
        v[i] = _f+i*_d;
      return  Sampling(v,_t);
    }
  }

  /**
   * Decimates this sampling. Beginning with the first sample, keeps only 
   * every m'th sample, while discarding the others in this sampling. If 
   * this sampling has n values, the new sampling will have 1+(n-1)/m values.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @param m the factor by which to decimate; must be positive.
   * @return the new sampling.
   */
  Sampling Sampling:: decimate(int m) {
    assert(m>0);
    int n = 1+(_n-1)/m;
    if (_v.size() == 0) {
      return  Sampling(n,m*_d,_f,_t);
    } else {
      vector<double> v(n);
      for (int i=0,j=0; i<n; ++i,j+=m)
        v[i] = _v[j];
      return Sampling(v,_t);
    }
  }

  /**
   * Interpolates this sampling. Inserts m-1 evenly spaced samples between 
   * each of the samples in this sampling. If this sampling has n values, 
   * the new sampling will have 1+(n-1)*m values.
   * <p>
   * This method returns a new sampling; it does not modify this sampling.
   * @param m the factor by which to interpolate.
   * @return the new sampling.
   */
  Sampling Sampling::interpolate(int m) {
    assert(m>0);
    int n = _n+(_n-1)*(m-1);
    if (_v.size() == 0) {
      return Sampling(n,_d/m,_f,_t);
    } else {
      vector<double> v (n);
      v[0] = _v[0];
      for (int i=1,j=m; i<_n; ++i,j+=m) {
        v[j] = _v[i];
        double dk = (v[j]-v[j-m])/m;
        double vk = v[j-m];
        for (int k=j-m+1; k<j; ++k,vk+=dk)
          v[k] = vk;
      }
      return Sampling(v,_t);
    }
  }

  double Sampling::value(int i) {
    return (_v.size() != 0)?_v[i]:_f+i*_d;
  }

  bool Sampling::almostEqual(double v1, double v2, double tiny) {
    double diff = v1-v2;
    return (diff<0.0)?-diff<tiny:diff<tiny;
  }

  double Sampling::tinyWith(Sampling s) {
    return (_td<s._td)?_td:s._td;
  }
