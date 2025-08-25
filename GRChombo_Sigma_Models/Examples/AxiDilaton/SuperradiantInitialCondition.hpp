#ifndef SUPERRADIANT_HPP_
#define SUPERRADIANT_HPP_

#include "Complex.hpp"
#include "SpheroidalHarmonics.hpp"

namespace SuperradiantInitCond
{
struct Superradiant
{
    double Real;
    double Im;
    double magnitude;
    Complex complex_value;
};

// auxiliary functions

inline Complex alpha_rad(const int n, const int em, const Complex freq,
                         const double spin)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex c0 = Complex(1.) - Complex(2.) * I * freq -
                 2. / b * I * (-spin * em / 2. + freq);

    return Complex(double(n * n)) + Complex(double(n)) * (Complex(1.) + c0) +
           c0;
}

inline Complex beta_rad(const int n, const int em, const Complex freq,
                        const double spin, const double mass,
                        const Complex sep_constant)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex q = -sqrt(mass * mass - freq * freq);

    Complex c1 = -Complex(4.) +
                 Complex(4.) * I * (freq - (b + Complex(1.)) * I * q) +
                 Complex(4.) / b * I * (-spin * em / 2. + freq) -
                 Complex(2.) * (freq * freq + q * q) / q;
    Complex c3 = -Complex(1.) +
                 Complex(2.) * I * powC(freq - I * q, Complex(3.)) / q +
                 b * powC(freq - I * q, Complex(2.)) * Complex(2.) +
                 spin * spin * q * q + 2. * spin * em * I * q - sep_constant -
                 (freq - I * q) * (freq - I * q) / q + b * q * Complex(2.) +
                 Complex(2.) / b * I *
                     (Complex(1.) + (freq - I * q) * (freq - I * q) / q) *
                     (-spin * em / 2. + freq);

    return -Complex(2. * n * n) + Complex(double(n)) * (Complex(2.) + c1) + c3;
}

inline Complex gamma_rad(const int n, const int em, const Complex freq,
                         const double spin, const double mass,
                         const Complex sep_constant)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex q = -sqrt(mass * mass - freq * freq);

    Complex c2 = Complex(3.) - Complex(2.) * I * freq -
                 Complex(2.) * (q * q - freq * freq) / q -
                 Complex(2.) / b * I * (-spin * em / 2. + freq);
    Complex c4 = powC(freq - I * q, Complex(4.)) / (q * q) +
                 Complex(2.) * I * freq * (freq - I * q) * (freq - I * q) / q -
                 Complex(2.) / b * I * (freq - I * q) * (freq - I * q) / q *
                     (-spin * em / 2. + freq);

    return Complex(double(n * n)) + Complex(double(n)) * (-Complex(3.) + c2) +
           c4;
}

inline Complex An_rad(const int n, const int em, const Complex freq,
                      const double spin, const double mass,
                      const Complex sep_constant)
{
    if (n == 0)
    {
        Complex result(1., 0);
        return result;
    }
    if (n == 1)
    {
        return -beta_rad(0, em, freq, spin, mass, sep_constant) /
               alpha_rad(0, em, freq, spin);
    }
    else
    {
        return (-beta_rad(n - 1, em, freq, spin, mass, sep_constant) *
                    An_rad(n - 1, em, freq, spin, mass, sep_constant) -
                gamma_rad(n - 1, em, freq, spin, mass, sep_constant) *
                    An_rad(n - 2, em, freq, spin, mass, sep_constant)) /
               alpha_rad(n - 1, em, freq, spin);
    }
}

inline Complex Radial_Component(const int em, const Complex freq,
                                const double spin, const double kerr_mass,
                                const double scalar_mass,
                                const Complex sep_constant, const double rr)
{
    Complex I(0., 1.);
    Complex q = -sqrt(scalar_mass * scalar_mass - freq * freq);
    double rplus = (kerr_mass + sqrt(kerr_mass * kerr_mass - spin * spin));
    double rminus = (kerr_mass - sqrt(kerr_mass * kerr_mass - spin * spin));
    double wcrit = em * spin / (2. * rplus);
    Complex sigma =
        Complex(1. / (rplus - rminus) * rplus * 2.) * (-Complex(wcrit) + freq);
    Complex chi_init =
        (scalar_mass * scalar_mass - Complex(2.) * freq * freq) / q;

    double rrBL = rr * (1. + rplus / (4. * rr)) * (1. + rplus / (4. * rr));

    Complex prefactor(0.0);
    if (rrBL > rplus)
    {
        prefactor =
            powC(Complex(rrBL - rplus), -I * sigma) *
            powC(Complex(rrBL - rminus), -Complex(1.) + I * sigma + chi_init) *
            exp(rrBL * q);
    }

    Complex sum;

    for (int i = 0; i <= 15; i++)
    {
        if (rrBL > rplus)
        {
            sum =
                sum + pow((rrBL - rplus) / (rrBL - rminus), i) *
                          An_rad(i, em, freq, spin, scalar_mass, sep_constant);
        }
    }

    return prefactor * sum;
}

inline Superradiant InitialConditionFunction(const double x, const double y,
                                             const double z,
                                             const double &scalarM,
                                             const double &kerrA,
                                             const double &kerrM)
{
    Superradiant result;

    Complex I(0., 1.);

    int em = 1;
    Complex freq(0.286675, 1.54269e-8);
    // double sep_const_Re = 2.00031;
    // double sep_const_Im = 1.43279e-9;
    Complex sep_constant(2.00031, 1.43279e-9);
    double Norm = 0.865918;

    double rr = simd_max(sqrt(x * x + y * y + z * z), 1e-6);

    SpheroidalHarmonics::S_lm_t S_lm011 = SpheroidalHarmonics::scalar_S_lm(
        x, y, z, em, freq, kerrA, scalarM, sep_constant, Norm);

    double phi_angle = atan2(y, x);

    // double kerra = kerrA;
    // double kerrm = kerrM;
    // double scalarm =scalarM;

    Complex RR =
        Radial_Component(em, freq, kerrA, kerrM, scalarM, sep_constant, rr);

    Complex temp = exp(-I * em * phi_angle) * RR * S_lm011.complex_value;

    result.Real = temp.getReal();     // /5000.;
    result.Im = temp.getImaginary();  // /5000.;
    result.magnitude = temp.getAbs(); // /5000.;
    result.complex_value = temp;      // /5000.;

    return result;
}

} // namespace SuperradiantInitCond

#endif /* SUPERRADIANT_HPP_ */

// namespace SuperradiantInitCond
// {
// template<class data_t> struct Superradiant
// {
//   data_t Real;
//   data_t Im;
//   data_t magnitude;
//   Complex<data_t> complex_value;
// };

// // auxiliary functions

// template <class data_t>
// Complex<data_t> alpha_rad(const int n, const int em, const Complex<data_t>
// freq, const data_t spin)
// {
//   Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1.  );
//   Complex<data_t> c0 = Complex<data_t>(1.) - Complex<data_t>(2.) * I * freq
//   - 2./ b * I  * ( - spin * em / 2.+freq );

//   return Complex<data_t>(double(n * n)) + Complex<data_t>(double(n)) *
//   (Complex<data_t>(1.)+ c0 )  + c0;
// }

// template <class data_t>
// Complex<data_t> beta_rad(const int n, const int em, const Complex<data_t>
// freq,
//               const data_t spin, const data_t mass, const Complex<data_t>
//               sep_constant)
// {
//   Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1. );
//   Complex<data_t> q = - sqrt(mass * mass - freq * freq);

//   Complex<data_t> c1 = - Complex<data_t>(4.) + Complex<data_t>(4.) * I *
//   (freq -  (b + Complex<data_t>(1.)) *I * q )
//                         + Complex<data_t>(4.)/ b  *I* (- spin * em/2. +freq)
//                         - Complex<data_t>(2.) *(freq*freq + q*q)/q;
//   Complex<data_t> c3 = - Complex<data_t>(1.)+ Complex<data_t>(2.)* I
//   *powC(freq - I *q, Complex<data_t>(3.))/q
//                             + b * powC(freq - I * q,
//                             Complex<data_t>(2.))*Complex<data_t>(2.) + spin *
//                             spin * q * q
//                             + 2. * spin * em * I * q   - sep_constant - (freq
//                             - I * q)*(freq - I * q)/q + b * q*
//                             Complex<data_t>(2.)  +
//                            Complex<data_t>(2.) / b * I  *
//                            (Complex<data_t>(1.)+(freq - I * q)*(freq - I *
//                            q)/q ) * (- spin * em/2.+freq );

//   return -Complex<data_t>(2. * n * n) + Complex<data_t>(double(n)) *
//   (Complex<data_t>(2.)+c1) + c3;
// }

// template <class data_t>
// Complex<data_t> gamma_rad(const int n, const int em, const Complex<data_t>
// freq,
//   const data_t spin, const data_t mass, const Complex<data_t> sep_constant)
// {
//     Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1.  );
//   Complex<data_t> q = - sqrt(mass * mass - freq * freq);

//   Complex<data_t> c2 = Complex<data_t>(3.) - Complex<data_t>(2.) * I * freq
//                       - Complex<data_t>(2.) * (q * q - freq * freq)/q -
//                       Complex<data_t>(2.)/ b  *I* (- spin * em/2.+freq );
//   Complex<data_t> c4 = powC(freq - I * q, Complex<data_t>(4.))/(q * q) +
//   Complex<data_t>(2.) * I * freq * (freq - I * q) * (freq - I * q)/q -
//           Complex<data_t>(2.) / b * I * (freq - I * q) * (freq - I * q) /q *
//           (- spin * em/2.+freq);

//   return Complex<data_t>(double(n * n)) + Complex<data_t>(double(n))*(-
//   Complex<data_t>(3.)+c2 ) + c4;
// }

// template <class data_t>
// Complex<data_t> An_rad(const int n, const int em, const Complex<data_t> freq,
// const data_t spin, const data_t mass, const Complex<data_t> sep_constant)
// {
//     if(n==0)
//         {Complex<data_t> result(1.,0);
//           return result;}
//     if(n==1)
//         {return - beta_rad(0, em, freq, spin, mass,
//         sep_constant)/alpha_rad(0, em, freq, spin); }
//     else
//         {return (-beta_rad(n-1, em, freq, spin, mass,
//         sep_constant)*An_rad(n-1, em, freq, spin, mass, sep_constant)
//                 - gamma_rad(n-1, em, freq, spin, mass,
//                 sep_constant)*An_rad(n-2, em, freq, spin, mass,
//                 sep_constant))/alpha_rad(n-1, em, freq, spin);
//         }
// }

// template <class data_t>
// Complex<data_t> Radial_Component(const int em, const Complex<data_t> freq,
// const data_t spin, const data_t kerr_mass, const data_t scalar_mass, const
// Complex<data_t> sep_constant, const data_t rr)
// {
//   Complex<data_t> I(0.,1.);
//   Complex<data_t> q = - sqrt(scalar_mass * scalar_mass - freq * freq);
//   data_t rplus = kerr_mass + sqrt( kerr_mass *  kerr_mass - spin * spin);
//   data_t rminus = kerr_mass - sqrt( kerr_mass * kerr_mass - spin * spin);
//   data_t wcrit = em * spin /(2. * rplus) ;
//   Complex<data_t> sigma = Complex<data_t>(1./(rplus - rminus)*rplus * 2.)* (-
//   Complex<data_t>(wcrit) + freq); Complex<data_t> chi_init = (scalar_mass *
//   scalar_mass - Complex<data_t>(2.) * freq * freq)/q;

//   Complex<data_t> prefactor = powC(Complex<data_t>(rr - rplus),-I * sigma)*
//                               powC(Complex<data_t>(rr - rminus), -
//                               Complex<data_t>(1.) + I* sigma + chi_init ) *
//                               exp(rr * q);

//   Complex<data_t> sum;

//   for (int i = 0; i <= 10; i++)
//   {
//       sum = sum + pow((rr - rplus)/(rr - rminus), i ) * An_rad(i, em, freq,
//       spin, scalar_mass, sep_constant);
//   }

//   return prefactor * sum;
// }

// template <class data_t>
// Superradiant<data_t> InitialConditionFunction(const data_t x, const double y,
// const double z, const data_t& scalarM, const data_t& kerrA, const data_t&
// kerrM)
// {
//     Superradiant<data_t> result;

//     int em = 1;
//     Complex<data_t> freq(0.286675, 1.54269e-8) ;
//     // data_t sep_const_Re = 2.00031;
//     // data_t sep_const_Im = 1.43279e-9;
//     Complex<data_t> sep_constant(2.00031, 1.43279e-9);
//     data_t Norm = 0.865918;

//     data_t rr = simd_max(sqrt(x * x + y * y + z * z), 1e-6);

//     SpheroidalHarmonics::S_lm_t<data_t>  S_lm011 =
//         SpheroidalHarmonics::scalar_S_lm(x, y, z, em, freq, kerrA, scalarM,
//         sep_constant, Norm);

//     data_t phi_angle = atan2(y, x);

//     // data_t kerra = kerrA;
//     // data_t kerrm = kerrM;
//     // data_t scalarm =scalarM;

//     Complex<data_t> RR = Radial_Component(em, freq, kerrA, kerrM, scalarM,
//     sep_constant, rr);

//     Complex<data_t> temp = exp(-em*phi_angle)* RR * S_lm011.complex_value;

//     result.Real = temp.getReal();
//     result.Im = temp.getImaginary();
//     result.magnitude = temp.getAbs();
//     result.complex_value = temp;

//     return result;
// }

// } // namespace SuperradiantInitCond

// #endif /* SUPERRADIANT_HPP_ */

#ifndef SUPERRADIANT_HPP_
#define SUPERRADIANT_HPP_

#include "Complex.hpp"
#include "SpheroidalHarmonics.hpp"

namespace SuperradiantInitCond
{
struct Superradiant
{
    double Real;
    double Im;
    double magnitude;
    Complex complex_value;
};

// auxiliary functions

inline Complex alpha_rad(const int n, const int em, const Complex freq,
                         const double spin)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex c0 = Complex(1.) - Complex(2.) * I * freq -
                 2. / b * I * (-spin * em / 2. + freq);

    return Complex(double(n * n)) + Complex(double(n)) * (Complex(1.) + c0) +
           c0;
}

inline Complex beta_rad(const int n, const int em, const Complex freq,
                        const double spin, const double mass,
                        const Complex sep_constant)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex q = -sqrt(mass * mass - freq * freq);

    Complex c1 = -Complex(4.) +
                 Complex(4.) * I * (freq - (b + Complex(1.)) * I * q) +
                 Complex(4.) / b * I * (-spin * em / 2. + freq) -
                 Complex(2.) * (freq * freq + q * q) / q;
    Complex c3 = -Complex(1.) +
                 Complex(2.) * I * powC(freq - I * q, Complex(3.)) / q +
                 b * powC(freq - I * q, Complex(2.)) * Complex(2.) +
                 spin * spin * q * q + 2. * spin * em * I * q - sep_constant -
                 (freq - I * q) * (freq - I * q) / q + b * q * Complex(2.) +
                 Complex(2.) / b * I *
                     (Complex(1.) + (freq - I * q) * (freq - I * q) / q) *
                     (-spin * em / 2. + freq);

    return -Complex(2. * n * n) + Complex(double(n)) * (Complex(2.) + c1) + c3;
}

inline Complex gamma_rad(const int n, const int em, const Complex freq,
                         const double spin, const double mass,
                         const Complex sep_constant)
{
    Complex I(0., 1.);
    double b = sqrt(-spin * spin + 1.);
    Complex q = -sqrt(mass * mass - freq * freq);

    Complex c2 = Complex(3.) - Complex(2.) * I * freq -
                 Complex(2.) * (q * q - freq * freq) / q -
                 Complex(2.) / b * I * (-spin * em / 2. + freq);
    Complex c4 = powC(freq - I * q, Complex(4.)) / (q * q) +
                 Complex(2.) * I * freq * (freq - I * q) * (freq - I * q) / q -
                 Complex(2.) / b * I * (freq - I * q) * (freq - I * q) / q *
                     (-spin * em / 2. + freq);

    return Complex(double(n * n)) + Complex(double(n)) * (-Complex(3.) + c2) +
           c4;
}

inline Complex An_rad(const int n, const int em, const Complex freq,
                      const double spin, const double mass,
                      const Complex sep_constant)
{
    if (n == 0)
    {
        Complex result(1., 0);
        return result;
    }
    if (n == 1)
    {
        return -beta_rad(0, em, freq, spin, mass, sep_constant) /
               alpha_rad(0, em, freq, spin);
    }
    else
    {
        return (-beta_rad(n - 1, em, freq, spin, mass, sep_constant) *
                    An_rad(n - 1, em, freq, spin, mass, sep_constant) -
                gamma_rad(n - 1, em, freq, spin, mass, sep_constant) *
                    An_rad(n - 2, em, freq, spin, mass, sep_constant)) /
               alpha_rad(n - 1, em, freq, spin);
    }
}

inline Complex Radial_Component(const int em, const Complex freq,
                                const double spin, const double kerr_mass,
                                const double scalar_mass,
                                const Complex sep_constant, const double rr)
{
    Complex I(0., 1.);
    Complex q = -sqrt(scalar_mass * scalar_mass - freq * freq);
    double rplus = (kerr_mass + sqrt(kerr_mass * kerr_mass - spin * spin));
    double rminus = (kerr_mass - sqrt(kerr_mass * kerr_mass - spin * spin));
    double wcrit = em * spin / (2. * rplus);
    Complex sigma =
        Complex(1. / (rplus - rminus) * rplus * 2.) * (-Complex(wcrit) + freq);
    Complex chi_init =
        (scalar_mass * scalar_mass - Complex(2.) * freq * freq) / q;

    double rrBL = rr * (1. + rplus / (4. * rr)) * (1. + rplus / (4. * rr));

    Complex prefactor(0.0);
    if (rrBL > rplus)
    {
        prefactor =
            powC(Complex(rrBL - rplus), -I * sigma) *
            powC(Complex(rrBL - rminus), -Complex(1.) + I * sigma + chi_init) *
            exp(rrBL * q);
    }

    Complex sum;

    for (int i = 0; i <= 15; i++)
    {
        if (rrBL > rplus)
        {
            sum =
                sum + pow((rrBL - rplus) / (rrBL - rminus), i) *
                          An_rad(i, em, freq, spin, scalar_mass, sep_constant);
        }
    }

    return prefactor * sum;
}

inline Superradiant InitialConditionFunction(const double x, const double y,
                                             const double z,
                                             const double &scalarM,
                                             const double &kerrA,
                                             const double &kerrM)
{
    Superradiant result;

    Complex I(0., 1.);

    int em = 1;
    Complex freq(0.286675, 1.54269e-8);
    // double sep_const_Re = 2.00031;
    // double sep_const_Im = 1.43279e-9;
    Complex sep_constant(2.00031, 1.43279e-9);
    double Norm = 0.865918;

    double rr = simd_max(sqrt(x * x + y * y + z * z), 1e-6);

    SpheroidalHarmonics::S_lm_t S_lm011 = SpheroidalHarmonics::scalar_S_lm(
        x, y, z, em, freq, kerrA, scalarM, sep_constant, Norm);

    double phi_angle = atan2(y, x);

    // double kerra = kerrA;
    // double kerrm = kerrM;
    // double scalarm =scalarM;

    Complex RR =
        Radial_Component(em, freq, kerrA, kerrM, scalarM, sep_constant, rr);

    Complex temp = exp(-I * em * phi_angle) * RR * S_lm011.complex_value;

    result.Real = temp.getReal() / 5000.;
    result.Im = temp.getImaginary() / 5000.;
    result.magnitude = temp.getAbs() / 5000.;
    result.complex_value = temp / 5000.;

    return result;
}

} // namespace SuperradiantInitCond

#endif /* SUPERRADIANT_HPP_ */

// namespace SuperradiantInitCond
// {
// template<class data_t> struct Superradiant
// {
//   data_t Real;
//   data_t Im;
//   data_t magnitude;
//   Complex<data_t> complex_value;
// };

// // auxiliary functions

// template <class data_t>
// Complex<data_t> alpha_rad(const int n, const int em, const Complex<data_t>
// freq, const data_t spin)
// {
//   Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1.  );
//   Complex<data_t> c0 = Complex<data_t>(1.) - Complex<data_t>(2.) * I * freq
//   - 2./ b * I  * ( - spin * em / 2.+freq );

//   return Complex<data_t>(double(n * n)) + Complex<data_t>(double(n)) *
//   (Complex<data_t>(1.)+ c0 )  + c0;
// }

// template <class data_t>
// Complex<data_t> beta_rad(const int n, const int em, const Complex<data_t>
// freq,
//               const data_t spin, const data_t mass, const Complex<data_t>
//               sep_constant)
// {
//   Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1. );
//   Complex<data_t> q = - sqrt(mass * mass - freq * freq);

//   Complex<data_t> c1 = - Complex<data_t>(4.) + Complex<data_t>(4.) * I *
//   (freq -  (b + Complex<data_t>(1.)) *I * q )
//                         + Complex<data_t>(4.)/ b  *I* (- spin * em/2. +freq)
//                         - Complex<data_t>(2.) *(freq*freq + q*q)/q;
//   Complex<data_t> c3 = - Complex<data_t>(1.)+ Complex<data_t>(2.)* I
//   *powC(freq - I *q, Complex<data_t>(3.))/q
//                             + b * powC(freq - I * q,
//                             Complex<data_t>(2.))*Complex<data_t>(2.) + spin *
//                             spin * q * q
//                             + 2. * spin * em * I * q   - sep_constant - (freq
//                             - I * q)*(freq - I * q)/q + b * q*
//                             Complex<data_t>(2.)  +
//                            Complex<data_t>(2.) / b * I  *
//                            (Complex<data_t>(1.)+(freq - I * q)*(freq - I *
//                            q)/q ) * (- spin * em/2.+freq );

//   return -Complex<data_t>(2. * n * n) + Complex<data_t>(double(n)) *
//   (Complex<data_t>(2.)+c1) + c3;
// }

// template <class data_t>
// Complex<data_t> gamma_rad(const int n, const int em, const Complex<data_t>
// freq,
//   const data_t spin, const data_t mass, const Complex<data_t> sep_constant)
// {
//     Complex<data_t> I(0.,1.);
//   data_t b = sqrt(- spin * spin+1.  );
//   Complex<data_t> q = - sqrt(mass * mass - freq * freq);

//   Complex<data_t> c2 = Complex<data_t>(3.) - Complex<data_t>(2.) * I * freq
//                       - Complex<data_t>(2.) * (q * q - freq * freq)/q -
//                       Complex<data_t>(2.)/ b  *I* (- spin * em/2.+freq );
//   Complex<data_t> c4 = powC(freq - I * q, Complex<data_t>(4.))/(q * q) +
//   Complex<data_t>(2.) * I * freq * (freq - I * q) * (freq - I * q)/q -
//           Complex<data_t>(2.) / b * I * (freq - I * q) * (freq - I * q) /q *
//           (- spin * em/2.+freq);

//   return Complex<data_t>(double(n * n)) + Complex<data_t>(double(n))*(-
//   Complex<data_t>(3.)+c2 ) + c4;
// }

// template <class data_t>
// Complex<data_t> An_rad(const int n, const int em, const Complex<data_t> freq,
// const data_t spin, const data_t mass, const Complex<data_t> sep_constant)
// {
//     if(n==0)
//         {Complex<data_t> result(1.,0);
//           return result;}
//     if(n==1)
//         {return - beta_rad(0, em, freq, spin, mass,
//         sep_constant)/alpha_rad(0, em, freq, spin); }
//     else
//         {return (-beta_rad(n-1, em, freq, spin, mass,
//         sep_constant)*An_rad(n-1, em, freq, spin, mass, sep_constant)
//                 - gamma_rad(n-1, em, freq, spin, mass,
//                 sep_constant)*An_rad(n-2, em, freq, spin, mass,
//                 sep_constant))/alpha_rad(n-1, em, freq, spin);
//         }
// }

// template <class data_t>
// Complex<data_t> Radial_Component(const int em, const Complex<data_t> freq,
// const data_t spin, const data_t kerr_mass, const data_t scalar_mass, const
// Complex<data_t> sep_constant, const data_t rr)
// {
//   Complex<data_t> I(0.,1.);
//   Complex<data_t> q = - sqrt(scalar_mass * scalar_mass - freq * freq);
//   data_t rplus = kerr_mass + sqrt( kerr_mass *  kerr_mass - spin * spin);
//   data_t rminus = kerr_mass - sqrt( kerr_mass * kerr_mass - spin * spin);
//   data_t wcrit = em * spin /(2. * rplus) ;
//   Complex<data_t> sigma = Complex<data_t>(1./(rplus - rminus)*rplus * 2.)* (-
//   Complex<data_t>(wcrit) + freq); Complex<data_t> chi_init = (scalar_mass *
//   scalar_mass - Complex<data_t>(2.) * freq * freq)/q;

//   Complex<data_t> prefactor = powC(Complex<data_t>(rr - rplus),-I * sigma)*
//                               powC(Complex<data_t>(rr - rminus), -
//                               Complex<data_t>(1.) + I* sigma + chi_init ) *
//                               exp(rr * q);

//   Complex<data_t> sum;

//   for (int i = 0; i <= 10; i++)
//   {
//       sum = sum + pow((rr - rplus)/(rr - rminus), i ) * An_rad(i, em, freq,
//       spin, scalar_mass, sep_constant);
//   }

//   return prefactor * sum;
// }

// template <class data_t>
// Superradiant<data_t> InitialConditionFunction(const data_t x, const double y,
// const double z, const data_t& scalarM, const data_t& kerrA, const data_t&
// kerrM)
// {
//     Superradiant<data_t> result;

//     int em = 1;
//     Complex<data_t> freq(0.286675, 1.54269e-8) ;
//     // data_t sep_const_Re = 2.00031;
//     // data_t sep_const_Im = 1.43279e-9;
//     Complex<data_t> sep_constant(2.00031, 1.43279e-9);
//     data_t Norm = 0.865918;

//     data_t rr = simd_max(sqrt(x * x + y * y + z * z), 1e-6);

//     SpheroidalHarmonics::S_lm_t<data_t>  S_lm011 =
//         SpheroidalHarmonics::scalar_S_lm(x, y, z, em, freq, kerrA, scalarM,
//         sep_constant, Norm);

//     data_t phi_angle = atan2(y, x);

//     // data_t kerra = kerrA;
//     // data_t kerrm = kerrM;
//     // data_t scalarm =scalarM;

//     Complex<data_t> RR = Radial_Component(em, freq, kerrA, kerrM, scalarM,
//     sep_constant, rr);

//     Complex<data_t> temp = exp(-em*phi_angle)* RR * S_lm011.complex_value;

//     result.Real = temp.getReal();
//     result.Im = temp.getImaginary();
//     result.magnitude = temp.getAbs();
//     result.complex_value = temp;

//     return result;
// }

// } // namespace SuperradiantInitCond

// #endif /* SUPERRADIANT_HPP_ */
