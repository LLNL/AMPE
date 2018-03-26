/***
 *  For the complex step derivative method:
 *  f'(x) ~  Im [ f(x+ih) ] / h
 *  Define a double complex class that inherits from the
 *  library complex type and overloads appropriate operators.
 *  Mon Jan  8 22:42:20 PST 2001
 ***/

#ifndef HDRcomplexify
#define HDRcomplexify

#include <complex>
using namespace std;

#ifndef HDRderivify
inline double real(const double& r) {
  /***
   *  So the real() statement can be used even with
   *  the double version of the code to be complexified.
   *  Most useful inside printf statements.
   ***/
return r;
}

inline double imag(const double& r) {
  return 0.;
}
#endif // HDRderivify

class cplx : public complex<double> {
public:
  cplx() : complex<double>() {};
  cplx(const double& d) : complex<double>(d) {};
  cplx(const double& r, const double& i) : complex<double>(r,i) {};
  cplx(const complex<double>& z) : complex<double>(z) {};
  cplx(const complex<float>& z) : complex<double>(z) {};
  ~cplx(void) {/*cout<<"destroy cplx!";*/}
  operator double() {return this->real();}
  operator int() {return int(this->real());}
  // relational operators
  // Conversion constructor should be able to take care of the
  // operator== and != calls with double, but MIPS compiler 
  // complains of ambiguous inheritance.  This should be more
  // efficient anyway.  (A hint of what comes below.)
  friend inline bool operator==(const cplx&,const cplx&);
  friend inline bool operator==(const cplx&,const double&);
  friend inline bool operator==(const double&,const cplx&);
  friend inline bool operator!=(const cplx&,const cplx&);
  friend inline bool operator!=(const cplx&,const double&);
  friend inline bool operator!=(const double&,const cplx&);
  friend inline bool operator>(const cplx&,const cplx&);
  friend inline bool operator>(const cplx&,const double&);
  friend inline bool operator>(const double&,const cplx&);
  friend inline bool operator<(const cplx&,const cplx&);
  friend inline bool operator<(const cplx&,const double&);
  friend inline bool operator<(const double&,const cplx&);
  friend inline bool operator>=(const cplx&,const cplx&);
  friend inline bool operator>=(const cplx&,const double&);
  friend inline bool operator>=(const double&,const cplx&);
  friend inline bool operator<=(const cplx&,const cplx&);
  friend inline bool operator<=(const cplx&,const double&);
  friend inline bool operator<=(const double&,const cplx&);
  // here's the annoying thing:
  // Every function in class complex<double> that returns a
  // complex<double> causes ambiguities with function overloading
  // resolution because of the mix of types cplx and
  // complex<double> and double and int in math expressions.
  // So, although they are inherited, must redefine them
  // to return type cplx:
  // basic arithmetic
  inline cplx operator+() const;
  inline cplx operator+(const cplx&) const;
  inline cplx operator+(const double&) const;
  inline cplx operator+(const int&) const;
  inline friend cplx operator+(const double&, const cplx&);
  inline friend cplx operator+(const int&, const cplx&);
  inline cplx operator-() const;
  inline cplx operator-(const cplx&) const;
  inline cplx operator-(const double&) const;
  inline cplx operator-(const int&) const;
  inline friend cplx operator-(const double&, const cplx&);
  inline friend cplx operator-(const int&, const cplx&);
  inline cplx operator*(const cplx&) const;
  inline cplx operator*(const double&) const;
  inline cplx operator*(const int&) const;
  inline friend cplx operator*(const double&, const cplx&);
  inline friend cplx operator*(const int&, const cplx&);
  inline cplx operator/(const cplx&) const;
  inline cplx operator/(const double&) const;
  inline cplx operator/(const int&) const;
  inline friend cplx operator/(const double&, const cplx&);
  inline friend cplx operator/(const int&, const cplx&);

  inline cplx operator++(void);//added by Andrei Schaffer

  // from <math.h>
  inline friend cplx sin(const cplx&);
  inline friend cplx sinh(const cplx&);
  inline friend cplx cos(const cplx&);
  inline friend cplx cosh(const cplx&);
  inline friend cplx tan(const cplx&);
  inline friend cplx tanh(const cplx&);
  inline friend cplx log10(const cplx&);
  inline friend cplx log(const cplx&);
  inline friend cplx sqrt(const cplx&);
  inline friend cplx exp(const cplx&);
  inline friend cplx pow(const cplx&, const cplx&);
  inline friend cplx pow(const cplx&, const double&);
  inline friend cplx pow(const cplx&, const int&);
  inline friend cplx pow(const double&, const cplx&);
  inline friend cplx pow(const int&, const cplx&);
  // complex versions of these are not in standard library
  // or they need to be redefined:
  // (frexp, modf, and fmod have not been dealt with)
  inline friend cplx fabs(const cplx&);
  inline friend cplx asin(const cplx&);
  inline friend cplx acos(const cplx&);
  inline friend cplx atan(const cplx&);
  inline friend cplx atan2(const cplx&, const cplx&);
  inline friend cplx ceil(const cplx&);
  inline friend cplx floor(const cplx&);
  inline friend cplx ldexp(const cplx&, const int&);
};


inline bool operator==(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) == real(rhs);
}

inline bool operator==(const cplx& lhs, const double& rhs)
{
  return real(lhs) == rhs;
}

inline bool operator==(const double& lhs, const cplx& rhs)
{
  return lhs == real(rhs);
}

inline bool operator!=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) != real(rhs);
}

inline bool operator!=(const cplx& lhs, const double& rhs)
{
  return real(lhs) != rhs;
}

inline bool operator!=(const double& lhs, const cplx& rhs)
{
  return lhs != real(rhs);
}

inline bool operator>(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) > real(rhs);
}

inline bool operator>(const cplx& lhs, const double& rhs)
{
  return real(lhs) > rhs;
}

inline bool operator>(const double& lhs, const cplx& rhs)
{
  return lhs > real(rhs);
}

inline bool operator<(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) < real(rhs);
}

inline bool operator<(const cplx& lhs, const double& rhs)
{
  return real(lhs) < rhs;
}

inline bool operator<(const double& lhs, const cplx& rhs)
{
  return lhs < real(rhs);
}

inline bool operator>=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) >= real(rhs);
}

inline bool operator>=(const cplx& lhs, const double& rhs)
{
  return real(lhs) >= rhs;
}

inline bool operator>=(const double& lhs, const cplx& rhs)
{
  return lhs >= real(rhs);
}

inline bool operator<=(const cplx& lhs, const cplx& rhs)
{
  return real(lhs) <= real(rhs);
}

inline bool operator<=(const cplx& lhs, const double& rhs)
{
  return real(lhs) <= rhs;
}

inline bool operator<=(const double& lhs, const cplx& rhs)
{
  return lhs <= real(rhs);
}

inline cplx cplx::operator++(void)
{
  (*this)=(*this)+1;
  return complex<double>(*this);
}

inline cplx cplx::operator+() const 
{
  return +complex<double>(*this);
}

inline cplx cplx::operator+(const cplx& z) const
{
  return complex<double>(*this)+complex<double>(z);
}

inline cplx cplx::operator+(const double& r) const
{
  return complex<double>(*this)+r;
}

inline cplx cplx::operator+(const int& i) const
{
  return complex<double>(*this)+double(i);
}

inline cplx operator+(const double& r, const cplx& z)
{
  return r+complex<double>(z);
}

inline cplx operator+(const int& i, const cplx& z)
{
  return double(i)+complex<double>(z);
}

inline cplx cplx::operator-() const
{
  return -complex<double>(*this);
}

inline cplx cplx::operator-(const cplx& z) const
{
  return complex<double>(*this)-complex<double>(z);
}

inline cplx cplx::operator-(const double& r) const
{
  return complex<double>(*this)-r;
}

inline cplx cplx::operator-(const int& i) const
{
  return complex<double>(*this)-double(i);
}

inline cplx operator-(const double& r, const cplx& z)
{
  return r-complex<double>(z);
}

inline cplx operator-(const int& i, const cplx& z)
{
  return double(i)-complex<double>(z);
}

inline cplx cplx::operator*(const cplx& z) const
{
  return complex<double>(*this)*complex<double>(z);
}

inline cplx cplx::operator*(const double& r) const
{
  return complex<double>(*this)*r;
}

inline cplx cplx::operator*(const int& i) const
{
  return complex<double>(*this)*double(i);
}

inline cplx operator*(const double& r, const cplx& z)
{
  return r*complex<double>(z);
}

inline cplx operator*(const int& i, const cplx& z)
{
  return double(i)*complex<double>(z);
}

inline cplx cplx::operator/(const cplx& z) const
{
  return complex<double>(*this)/complex<double>(z);
}

inline cplx cplx::operator/(const double& r) const
{
  return complex<double>(*this)/r;
}

inline cplx cplx::operator/(const int& i) const
{
  return complex<double>(*this)/double(i);
}

inline cplx operator/(const double& r, const cplx& z)
{
  return r/complex<double>(z);
}

inline cplx operator/(const int& i, const cplx& z)
{
  return double(i)/complex<double>(z);
}

inline cplx sin(const cplx& z)
{
  return sin(complex<double>(z));
}

inline cplx sinh(const cplx& z)
{
  return sinh(complex<double>(z));
}

inline cplx cos(const cplx& z)
{
  return cos(complex<double>(z));
}

inline cplx cosh(const cplx& z)
{
  return cosh(complex<double>(z));
}

#ifdef __GNUC__ // bug in gcc ?? get segv w/egcs-2.91.66 and 2.95.2
inline cplx tan(const cplx& z) 
{
  return sin(complex<double>(z))/cos(complex<double>(z));
}

inline cplx tanh(const cplx& z)
{
  return sinh(complex<double>(z))/cosh(complex<double>(z));
}

inline cplx log10(const cplx& z)
{
  return log(complex<double>(z))/log(10.);
}
#else
inline cplx tan(const cplx& z)
{
  return tan(complex<double>(z));
}

inline cplx tanh(const cplx& z)
{
  return tanh(complex<double>(z));
}

inline cplx log10(const cplx& z)
{
  return log10(complex<double>(z));
}
#endif

inline cplx log(const cplx& z)
{
  return log(complex<double>(z));
}

inline cplx sqrt(const cplx& z)
{
  return sqrt(complex<double>(z));
}

inline cplx exp(const cplx& z)
{
  return exp(complex<double>(z));
}

inline cplx pow(const cplx& a, const cplx& b)
{
  return pow(complex<double>(a),complex<double>(b));
}

inline cplx pow(const cplx& a, const double& b)
{
  return pow(complex<double>(a),b);
}

inline cplx pow(const cplx& a, const int& b)
{
  return pow(complex<double>(a),double(b));
}

inline cplx pow(const double& a, const cplx& b)
{
  return pow(a,complex<double>(b));
}

inline cplx pow(const int& a, const cplx& b)
{
  return pow(double(a),complex<double>(b));
}

inline cplx fabs(const cplx& z)
{
  return (real(z)<0.0) ? -z:z;
}

#define surr_TEENY (1.e-24) /* machine zero compared to nominal magnitude of
			       the real part */

inline cplx asin(const cplx& z)
{
  // derivative trouble if imag(z) = +/- 1.0
  return cplx(asin(real(z)),imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

inline cplx acos(const cplx& z)
{
  // derivative trouble if imag(z) = +/- 1.0
  return cplx(acos(real(z)),-imag(z)/sqrt(1.0-real(z)*real(z)+surr_TEENY));
}

#undef surr_TEENY

inline cplx atan(const cplx& z)
{
  return cplx(atan(real(z)),imag(z)/(1.0+real(z)*real(z)));
}

inline cplx atan2(const cplx& z1, const cplx& z2)
{
  return cplx(atan2(real(z1),real(z2)),
	      (real(z2)*imag(z1)-real(z1)*imag(z2))
	      /(real(z1)*real(z1)+real(z2)*real(z2)));
}

inline cplx ceil(const cplx& z)
{
  return cplx(ceil(real(z)),0.);
}

inline cplx floor(const cplx& z)
{
  return cplx(floor(real(z)),0.);
}

inline cplx ldexp(const cplx& z, const int& i)
{
  return cplx(ldexp(real(z),i),ldexp(imag(z),i));
}


class cplx_vector
{
 public:
  void *re;
  void *im;

  cplx_vector(void *r, void *i): re(r),im(i){}

  cplx_vector(void): re(0),im(0){}

  void* real(){return re;}
  void* imag(){return im;}
};

#endif
