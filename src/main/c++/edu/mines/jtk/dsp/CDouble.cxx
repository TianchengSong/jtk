/*! \file CDouble.cxx */

/******************************************************************************
 * C++ port from Java by Tiancheng Song
 * Please refer to the header file for complete documentation and licensing information.
 *
 * Member Functions:
 * 
 * Cdouble::Cdouble()
 * Cdouble::Cdouble(double r)
 * Cdouble::Cdouble(double r, double i)
 * Cdouble::Cdouble(const Cdouble &x)
 * Cdouble::Cdouble plus(const Cdouble &x)
 * Cdouble::minus(Cdouble x);
 * Cdouble::times(Cdouble x);
 * Cdouble::over(Cdouble x);
 * Cdouble::plus(double x);
 * Cdouble::minus(double x);
 * Cdouble::times(double x);
 * Cdouble::over(double x);
 * Cdouble::plusEquals(Cdouble x);
 * Cdouble::minusEquals(Cdouble x);
 * Cdouble::timesEquals(Cdouble x);
 * Cdouble::overEquals(Cdouble x);
 * Cdouble::plusEquals(double x);
 * Cdouble::minusEquals(double x);
 * Cdouble::timesEquals(double x);
 * Cdouble::overEquals(double x);
 * Cdouble::conjEquals();
 * Cdouble::invEquals();
 * Cdouble::negEquals();
 * Cdouble::isReal();
 * Cdouble::isImag();
 * Cdouble::conj();
 * Cdouble::inv();
 * Cdouble::neg();
 * Cdouble::abs();
 * Cdouble::arg();
 * Cdouble::norm();
 * Cdouble::sqrt();
 * Cdouble::exp();
 * Cdouble::log();
 * Cdouble::log10();
 * Cdouble::pow(double y);
 * Cdouble::pow(Cdouble y);
 * Cdouble::sin();
 * Cdouble::cos();
 * Cdouble::tan();
 * Cdouble::sinh();
 * Cdouble::cosh();
 * Cdouble::tanh();
 * Cdouble::isReal(Cdouble x);
 * Cdouble::isImag(Cdouble x) ;
 * Cdouble::conj(Cdouble x);
 * Cdouble::inv(Cdouble x);
 * Cdouble::neg(Cdouble x);
 * Cdouble::polar(double r, double a);
 * Cdouble::add(Cdouble x, Cdouble y);
 * Cdouble::sub(Cdouble x, Cdouble y);
 * Cdouble::mul(Cdouble x, Cdouble y);
 * Cdouble::div(Cdouble x, Cdouble y);
 * Cdouble::abs(Cdouble x);
 * Cdouble::arg(Cdouble x);
 * Cdouble::norm(Cdouble x);
 * Cdouble::sqrt(Cdouble x);
 * Cdouble::exp(Cdouble x);
 * Cdouble::log(Cdouble x);
 * Cdouble::log10(Cdouble x);
 * Cdouble::pow(Cdouble x, double y);
 * Cdouble::pow(double x, Cdouble y);
 * Cdouble::pow(Cdouble x, Cdouble y);
 * Cdouble::sin(Cdouble x);
 * Cdouble::cos(Cdouble x);
 * Cdouble::tan(Cdouble x);
 * Cdouble::sinh(Cdouble x);
 * Cdouble::cosh(Cdouble x);
 * Cdouble::tanh(Cdouble x);
 * Cdouble::equals(Cdouble obj);
 * Cdouble::uint64_t doubleToRawBits(double x);
 * Cdouble::hashCode();
 * Cdouble::toString();
 * Cdouble::max(double x, double y);
 * Cdouble::abs(double x);
 * Cdouble::sqrt(double x);
 * Cdouble::cos(double x);
 * Cdouble::sin(double x);
 * Cdouble::cosh(double x);
 * Cdouble::sinh(double x);
 * Cdouble::exp(double x);
 * Cdouble::log(double x);
 * Cdouble::pow(double x, double y);
 * Cdouble::atan2(double y, double x);
 ******************************************************************************
*/

#include "CDouble.h"
#include <math.h>
using namespace std;
/**
  * Constructs a complex number with zero real and imaginary parts.
  */
Cdouble::Cdouble() {
  r = 0.0;
  i = 0.0;
}

/**
  * Constructs a complex number with zero imaginary part.
  * @param r the real part.
  */
Cdouble::Cdouble(double r) {
  this->r = r; 
  i = 0;;
}

/**
  * Constructs a complex number.
  * @param r the real part.
  * @param i the imaginary part.
  */
Cdouble::Cdouble(double r, double i) {
  this->r = r;
  this->i = i;
}

/**
  * Constructs a copy of the specified complex number.
  * @param x the complex number.
  */

Cdouble::Cdouble(const Cdouble &x) {
  this->r = x.r;
  this->i = x.i;    
}

/**
  * Returns the sum z + x, where z is this complex number.
  * @param x a complex number.
  * @return z + x.
  */
Cdouble Cdouble::plus(const Cdouble &x) {
  return Cdouble(*this).plusEquals(x);
}

/**
  * Returns the difference z - x, where z is this complex number.
  * @param x a complex number.
  * @return z - x.
  */
Cdouble Cdouble::minus(Cdouble x) {
  return Cdouble(*this).minusEquals(x);
}

/**
  * Returns the product z * x, where z is this complex number.
  * @param x a complex number.
  * @return z * x.
  */
Cdouble Cdouble::times(Cdouble x) {
  return Cdouble(*this).timesEquals(x);
}

/**
  * Returns the quotent z / x, where z is this complex number.
  * @param x a complex number.
  * @return z / x.
  */
Cdouble Cdouble::over(Cdouble x) {
  return Cdouble(*this).overEquals(x);
}

/**
  * Returns the sum z + x, where z is this complex number.
  * @param x a real number.
  * @return z + x.
  */
Cdouble Cdouble::plus(double x) {
  return Cdouble(*this).plusEquals(x);
}

/**
  * Returns the difference z - x, where z is this complex number.
  * @param x a real number.
  * @return z - x.
  */
Cdouble Cdouble::minus(double x) {
  return Cdouble(*this).minusEquals(x);
}

/**
  * Returns the product z * x, where z is this complex number.
  * @param x a real number.
  * @return z * x.
  */
Cdouble Cdouble::times(double x) {
  return Cdouble(*this).timesEquals(x);
}

/**
  * Returns the quotent z / x, where z is this complex number.
  * @param x a real number.
  * @return z / x.
  */
Cdouble Cdouble::over(double x) {
  return Cdouble(*this).overEquals(x);
}

/**
  * Returns the sum z += x, where z is this complex number.
  * @param x a complex number.
  * @return z += x.
  */
Cdouble Cdouble::plusEquals(Cdouble x) {
  r += x.r;
  i += x.i;
  return *this;
}

/**
  * Returns the difference z -= x, where z is this complex number.
  * @param x a complex number.
  * @return z -= x.
  */
Cdouble Cdouble::minusEquals(Cdouble x) {
  r -= x.r;
  i -= x.i;
  return *this;
}

/**
  * Returns the product z *= x, where z is this complex number.
  * @param x a complex number.
  * @return z *= x.
  */
Cdouble Cdouble::timesEquals(Cdouble x) {
  double tr = this->r;
  double ti = this->i;
  double xr = x.r;
  double xi = x.i;
  r = tr*xr-ti*xi;
  i = tr*xi+ti*xr;
  return *this;
}

/**
  * Returns the quotient z /= x, where z is this complex number.
  * @param x a complex number.
  * @return z /= x.
  */
Cdouble Cdouble::overEquals(Cdouble x) {
  double tr = this->r;
  double ti = this->i;
  double xr = x.r;
  double xi = x.i;
  double d = norm(x);
  r = (tr*xr+ti*xi)/d;
  i = (ti*xr-tr*xi)/d;
  return *this;
}

/**
  * Returns the sum z += x, where z is this complex number.
  * @param x a real number.
  * @return z += x.
  */
Cdouble Cdouble::plusEquals(double x) {
  r += x;
  return *this;
}

/**
  * Returns the difference z -= x, where z is this complex number.
  * @param x a real number.
  * @return z -= x.
  */
Cdouble Cdouble::minusEquals(double x) {
  r -= x;
  return *this;
}

/**
  * Returns the product z *= x, where z is this complex number.
  * @param x a real number.
  * @return z *= x.
  */
Cdouble Cdouble::timesEquals(double x) {
  r *= x;
  i *= x;
  return *this;
}

/**
  * Returns the quotient z /= x, where z is this complex number.
  * @param x a real number.
  * @return z /= x.
  */
Cdouble Cdouble::overEquals(double x) {
  r /= x;
  i /= x;
  return *this;
}

/**
  * Returns the conjugate z = conj(z), where z is this complex number.
  * @return z = conj(z).
  */
Cdouble Cdouble::conjEquals() {
  i = -i;
  return *this;
}

/**
  * Returns the inverse z = inv(z), where z is this complex number.
  * @return z = inv(z).
  */
Cdouble Cdouble::invEquals() {
  double d = norm();
  r =  r/d;
  i = -i/d;
  return *this;
}

/**
  * Returns the negative z = neg(z), where z is this complex number.
  * @return z = neg(z).
  */
Cdouble Cdouble::negEquals() {
  r = -r;
  i = -i;
  return *this;
}

/**
  * Determines whether this complex number is real (has zero imaginary part).
  * @return true, if real; false, otherwise.
  */
bool Cdouble::isReal() {
  return i==0.0;
}

/**
  * Determines whether this complex number is imaginary (has zero real part).
  * @return true, if imaginary; false, otherwise.
  */
bool Cdouble::isImag() {
  return r==0.0;
}

/**
  * Returns the complex conjugate of this complex number.
  * @return the complex conjugate.
  */
Cdouble Cdouble::conj() {
  return Cdouble(r,-i);
}

/**
  * Returns the complex inverse of this complex number.
  * @return the complex inverse.
  */
Cdouble Cdouble::inv() {
  double d = norm();
  return Cdouble(r/d,-i/d);
}

/**
  * Returns the complex negative of this complex number.
  * @return the complex negative.
  */
Cdouble Cdouble::neg() {
  return Cdouble(-r,-i);
}

/**
  * Returns the magnitude of this complex number.
  * @return the magnitude.
  */
double Cdouble::abs() {
  return abs(*this);
}

/**
  * Returns the argument of this complex number.
  * @return the argument.
  */
double Cdouble::arg() {
  return arg(*this);
}

/**
  * Returns the norm of this complex number.
  * The norm is the sum of the squares of the real and imaginary parts.
  * @return the norm.
  */
double Cdouble::norm() {
  return norm(*this);
}

/**
  * Returns the square-root of this complex number.
  * @return the square-root.
  */
Cdouble Cdouble::sqrt() {
  return sqrt(*this);
}

/**
  * Returns the exponential of this complex number.
  * @return the exponential.
  */
Cdouble Cdouble::exp() {
  return exp(*this);
}

/**
  * Returns the natural logarithm of this complex number.
  * @return the natural logarithm.
  */
Cdouble Cdouble::log() {
  return log(*this);
}

/**
  * Returns the logarithm base 10 of this complex number.
  * @return the logarithm base 10.
  */
Cdouble Cdouble::log10() {
  return log10(*this);
}

/**
  * Returns z to the y'th power, where z is this complex number.
  * @param y a real number.
  * @return z to the y'th power.
  */
Cdouble Cdouble::pow(double y) {
  return pow(*this,y);
}

/**
  * Returns z to the y'th power, where z is this complex number.
  * @param y a complex number.
  * @return z to the y'th power.
  */
Cdouble Cdouble::pow(Cdouble y) {
  return pow(*this,y);
}

/**
  * Returns the sine of this complex number.
  * @return the sine.
  */
Cdouble Cdouble::sin() {
  return sin(*this);
}

/**
  * Returns the cosine of this complex number.
  * @return the cosine.
  */
Cdouble Cdouble::cos() {
  return cos(*this);
}

/**
  * Returns the tangent of this complex number.
  * @return the tangent.
  */
Cdouble Cdouble::tan() {
  return tan(*this);
}

/**
  * Returns the hyberbolic sine of this complex number.
  * @return the hyberbolic sine.
  */
Cdouble Cdouble::sinh() {
  return sinh(*this);
}

/**
  * Returns the hyberbolic cosine of this complex number.
  * @return the hyberbolic cosine.
  */
Cdouble Cdouble::cosh() {
  return cosh(*this);
}

/**
  * Returns the hyberbolic tangent of this complex number.
  * @return the hyberbolic tangent.
  */
Cdouble Cdouble::tanh() {
  return tanh(*this);
}

/**
  * Determines whether x is real (has zero imaginary part).
  * @param x a complex number.
  * @return true, if real; false, otherwise.
  */
bool Cdouble::isReal(Cdouble x) {
  return x.i==0.0;
}

/**
  * Determines whether x is imaginary (has zero real part).
  * @param x a complex number.
  * @return true, if imaginary; false, otherwise.
  */
bool Cdouble::isImag(Cdouble x) {
  return x.r==0.0;
}

/**
  * Returns the conjugate of x.
  * @param x a complex number.
  * @return the conjugate.
  */
Cdouble Cdouble::conj(Cdouble x) {
  return Cdouble(x.r,-x.i);
}

/**
  * Returns the inverse of x.
  * @param x a complex number.
  * @return the complex inverse.
  */
Cdouble Cdouble::inv(Cdouble x) {
  double d = x.norm();
  return Cdouble(x.r/d,-x.i/d);
}

/**
  * Returns the negative of x.
  * @param x a complex number.
  * @return the negative.
  */
Cdouble Cdouble::neg(Cdouble x) {
  return Cdouble(-x.r,-x.i);
}

/**
  * Returns the complex number (r*cos(a),r*sin(a)).
  * @param r the polar radius.
  * @param a the polar angle.
  * @return the complex number.
  */
Cdouble Cdouble::polar(double r, double a) {
  return Cdouble(r*std::cos(a),r*std::sin(a));
}

/**
  * Returns the sum x + y.
  * @param x a complex number.
  * @param y a complex number.
  * @return the sum.
  */
Cdouble Cdouble::add(Cdouble x, Cdouble y) {
  return x.plus(y);
}

/**
  * Returns the difference x - y.
  * @param x a complex number.
  * @param y a complex number.
  * @return the difference.
  */
Cdouble Cdouble::sub(Cdouble x, Cdouble y) {
  return x.minus(y);
}

/**
  * Returns the product x * y.
  * @param x a complex number.
  * @param y a complex number.
  * @return the product.
  */
Cdouble Cdouble::mul(Cdouble x, Cdouble y) {
  return x.times(y);
}

/**
  * Returns the quotient x * y.
  * @param x a complex number.
  * @param y a complex number.
  * @return the quotient.
  */
Cdouble Cdouble::div(Cdouble x, Cdouble y) {
  return x.over(y);
}

/**
  * Returns the magnitude of a complex number.
  * @param x a complex number.
  * @return the magnitude.
  */
double Cdouble::abs(Cdouble x) {
  double ar = abs(x.r);
  double ai = abs(x.i);
  double s = max(abs(ar),abs(ai));
  if (s==0.0)
    return 0.0;
  ar /= s;
  ai /= s;
  return s*std::sqrt(ar*ar+ai*ai);
}

/**
  * Returns the argument of a complex number.
  * @param x a complex number.
  * @return the argument.
  */
double Cdouble::arg(Cdouble x) {
  return atan2(x.i,x.r);
}

/**
  * Returns the norm of a complex number.
  * The norm is the sum of the squares of the real and imaginary parts.
  * @param x a complex number.
  * @return the norm.
  */
double Cdouble::norm(Cdouble x) {
  return x.r*x.r+x.i*x.i;
}

/**
  * Returns the square root of a complex number.
  * @param x a complex number.
  * @return the square root.
  */
Cdouble Cdouble::sqrt(Cdouble x) {
  if (x.r==0.0) {
    double t = sqrt(0.5*abs(x.i));
    return Cdouble(t,(x.i<0.0)?-t:t);
  } else {
    double t = sqrt(2.0*(abs(x)+abs(x.r)));
    double u = 0.5*t;
    return (x.r>0.0) ? 
      Cdouble(u,x.i/t) :
      Cdouble(abs(x.i)/t,(x.i<0.0)?-u:u);
  }
}

/**
  * Returns the exponential of a complex number.
  * @param x a complex number.
  * @return the exponential.
  */
Cdouble Cdouble::exp(Cdouble x) {
  return polar(exp(x.r),x.i);
}

/**
  * Returns the natural logarithm of a complex number.
  * @param x a complex number.
  * @return the natural logarithm.
  */
Cdouble Cdouble::log(Cdouble x) {
  return Cdouble(log(abs(x)),arg(x));
}

/**
  * Returns the logarithm base 10 of a complex number.
  * @param x a complex number.
  * @return the logarithm base 10.
  */
Cdouble Cdouble::log10(Cdouble x) {
  return log(x).overEquals(log(10.0));
}

/**
  * Returns x to the y'th power.
  * @param x a complex number.
  * @param y a real number.
  * @return x to the y'th power.
  */
Cdouble Cdouble::pow(Cdouble x, double y) {
  if (x.i==0.0)
    return Cdouble(pow(x.r,y));
  Cdouble t = log(x);
  return polar(exp(y*t.r),y*t.i);
}

/**
  * Returns x to the y'th power.
  * @param x a real number.
  * @param y a complex number.
  * @return x to the y'th power.
  */
Cdouble Cdouble::pow(double x, Cdouble y) {
  if (x==0.0)
    return Cdouble();
  return polar(pow(x,y.r),y.i*log(x));
}

/**
  * Returns x to the y'th power.
  * @param x a complex number.
  * @param y a complex number.
  * @return x to the y'th power.
  */
Cdouble Cdouble::pow(Cdouble x, Cdouble y) {
  if (x.r==0.0 && x.i==0.0)
    return  Cdouble();
  return exp(y.times(log(x)));
}

/**
  * Returns the sine of a complex number.
  * @param x a complex number.
  * @return the sine.
  */
  Cdouble Cdouble::sin(Cdouble x) {
  return Cdouble(sin(x.r)*cosh(x.i),cos(x.r)*sinh(x.i));
}

/**
  * Returns the cosine of a complex number.
  * @param x a complex number.
  * @return the cosine.
  */
  Cdouble Cdouble::cos(Cdouble x) {
  return  Cdouble(cos(x.r)*cosh(x.i),-sin(x.r)*sinh(x.i));
}

/**
  * Returns the tangent of a complex number.
  * @param x a complex number.
  * @return the tangent.
  */
Cdouble Cdouble::tan(Cdouble x) {
  return sin(x).overEquals(cos(x));
}

/**
  * Returns the hyperbolic sine of a complex number.
  * @param x a complex number.
  * @return the hyperbolic sine.
  */
Cdouble Cdouble::sinh(Cdouble x) {
  return  Cdouble(sinh(x.r)*cos(x.i),cosh(x.r)*sin(x.i));
}

/**
  * Returns the hyperbolic cosine of a complex number.
  * @param x a complex number.
  * @return the hyperbolic cosine.
  */
Cdouble Cdouble::cosh(Cdouble x) {
  return  Cdouble(cosh(x.r)*cos(x.i),sinh(x.r)*sin(x.i));
}

/**
  * Returns the hyperbolic tangent of a complex number.
  * @param x a complex number.
  * @return the hyperbolic tangent.
  */
Cdouble Cdouble::tanh(Cdouble x) {
  return sinh(x).overEquals(cosh(x));
}

bool Cdouble::equals(Cdouble obj) {
  if (this == &obj)
    return true;
  if (&obj == NULL)
    return false;
  return (this->r==obj.r) && (this->i==obj.i);
}



inline uint64_t Cdouble::doubleToRawBits(double x) {
  uint64_t bits;
  memcpy(&bits, &x, sizeof bits);
  return bits;
}

int Cdouble::hashCode() {
  long rbits = doubleToRawBits(r);
  long ibits = doubleToRawBits(i);
  uint64_t rbitsU = (uint64_t)rbits;
  uint64_t ibitsU = (uint64_t)ibits;
  return (int)(rbitsU ^ (rbitsU >> 32) ^ ibitsU^ (ibitsU >> 32) );
}

string Cdouble::toString() {
  std::ostringstream s;
  if (i==0.0) {
    s << "(" << r << "+0.0i)";
  } else if (i>0.0) {
    s << "(" << r << "+" << i << "i)";
  } else {
    s << "(" << r << "-" << (-i) << "i)";
  }
  return s.str();
}

double Cdouble::max(double x, double y) {
  return (x>=y)?x:y;
}

double Cdouble::abs(double x) {
  return (x>=0.0)?x:-x;
}

double Cdouble::sqrt(double x) {
  return std::sqrt(x);  
}
double Cdouble::cos(double x) {
  return std::cos(x);  
}

double Cdouble::sin(double x) {
  return std::sin(x);  
}

double Cdouble::cosh(double x) {
  return std::cosh(x);  
}

double Cdouble::sinh(double x) {
  return std::sinh(x);  
}

double Cdouble::exp(double x) {
  return std::exp(x);  
}

double Cdouble::log(double x) {
  return std::log(x);
}

double Cdouble::pow(double x, double y) {
  return std::pow(x,y);
}

double Cdouble::atan2(double y, double x) {
  return std::atan2(y,x);
}


