// Interpolation function for 3 phases
// Folch and Plapp (2005)
// note: p1, p2 order does not matter
// Derivatives are constrained derivatives (p0+p1+p2=1)
double g0(const double p0, const double p1, const double p2)
{
   return 0.25 * p0 * p0 *
          (15. * (1. - p0) * (1. + p0 - (p2 - p1) * (p2 - p1)) +
           p0 * (9. * p0 * p0 - 5.));
}

// dg/dp0 - 1/3 * sum_j dg/dpj
// note: p1, p2 order does not matter
double dg0dp0(const double p0, const double p1, const double p2)
{
   return 2.5 * p0 *
          ((p2 - p1) * (p2 - p1) * (3. * p0 - 2.) +
           (1. - p0) * (1. - p0) * (3. * p0 + 2.));
}

double dg0dp1(const double p0, const double p1, const double p2)
{
   return -0.5 * dg0dp0(p1, p0, p2) + 7.5 * (p1 * p1 * (1. - p1) * (p2 - p0));
}

// Tilting function modification to stabilize Folch and Plapp function
// outside [0.,1.]
// Schwen, Jiang, Aagesen, Comput. Mat. Sci. 195 (2021), 110466
double soft_heaviside(const double phi)
{
   if (phi < 0.)
      return 0.;
   else if (phi > 1.)
      return 1.;

   return phi * phi * (3. - 2. * phi);
}

double dsoft_heaviside(const double phi)
{
   if (phi < 0.)
      return 0.;
   else if (phi > 1.)
      return 0.;

   return phi * 6. * (1. - phi);
}

double g0_soft_heaviside(const double p0, const double p1, const double p2)
{
   return g0(soft_heaviside(p0), soft_heaviside(p1), soft_heaviside(p2));
}

double dg0dp0_soft_heaviside(const double p0, const double p1, const double p2)
{
   //   return dg0dp0(soft_heaviside(p0), soft_heaviside(p1),
   //   soft_heaviside(p2)) *
   //          dsoft_heaviside(p0);
   double eps = 1.e-8;
   double val0 = g0_soft_heaviside(p0, p1, p2);
   double val1 =
       g0_soft_heaviside(p0 + 2. * eps / 3., p1 - eps / 3., p2 - eps / 3.);
   return (val1 - val0) / eps;
}

double dg0dp1_soft_heaviside(const double p0, const double p1, const double p2)
{
   //   return dg0dp1(soft_heaviside(p0), soft_heaviside(p1),
   //   soft_heaviside(p2))
   //        * dsoft_heaviside(p1);
   double eps = 1.e-8;
   double val0 = g0_soft_heaviside(p0, p1, p2);
   double val1 =
       g0_soft_heaviside(p0 - eps / 3., p1 + 2. * eps / 3., p2 - eps / 3.);
   return (val1 - val0) / eps;
}
