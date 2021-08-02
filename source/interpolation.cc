
// Interpolation function for 3 phases
// Folch and Plapp (2005)
double g0(const double p0, const double p1, const double p2)
{
   return 0.25 * p0 * p0 *
          (15. * (1. - p0) * (1. + p0 - (p2 - p1) * (p2 - p1)) +
           p0 * (9. * p0 * p0 - 5.));
}

// dg/dp0 - 1/3 * sum_j dg/dpj
double dg0dp0(const double p0, const double p1, const double p2)
{
   return 2.5 * p0 *
          ((p2 - p1) * (p2 - p1) * (3. * p0 - 2.) +
           (1. - p0) * (1. - p0) * (3. * p0 + 2.));
}

double dg0dp1(const double p0, const double p1, const double p2)
{
   return -0.5 * dg0dp0(p0, p1, p2) + 7.5 * (p0 * p0 * (1. - p0) * (p2 - p1));
}
