// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE. 
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
/*
 * Description: High-level wrapper for quaternion system solver
 */
#define SAMRAI_INLINE_KEYWORD inline

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setVerbose(const bool verbose)
{
  d_verbose = verbose;
  d_fac_ops->setVerbose(verbose);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setProlongationMethod(const std::string & prolongation_method)
{
  d_fac_ops->setProlongationMethod(prolongation_method);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setLevelSolverTolerance(double tol)
{
  d_fac_ops->setLevelSolverTolerance(tol);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setLevelSolverMaxIterations(int max_iterations)
{
  d_fac_ops->setLevelSolverMaxIterations(max_iterations);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setCoarsestLevelSolverTolerance(double tol)
{
  d_fac_ops->setCoarsestLevelSolverTolerance(tol);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setCoarsestLevelSolverMaxIterations(int max_iterations)
{
  d_fac_ops->setCoarsestLevelSolverMaxIterations(max_iterations);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setCoarseFineDiscretization(const std::string & coarsefine_method)
{
  d_fac_ops->setCoarseFineDiscretization(coarsefine_method);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setMaxCycles(int max_cycles) {
  d_fac_solver.setMaxCycles(max_cycles);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setResidualTolerance(double residual_tol) {
  d_fac_solver.setResidualTolerance(residual_tol);
}

SAMRAI_INLINE_KEYWORD int
QuatSysSolver::getNumberOfIterations() const
{
  return d_fac_solver.getNumberOfIterations();
}

SAMRAI_INLINE_KEYWORD double
QuatSysSolver::getResidualNorm() const
{
  return d_fac_solver.getResidualNorm();
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::getFACConvergenceFactors(double &   avg_factor,
					     double & final_factor) const
{
  d_fac_solver.getConvergenceFactors(avg_factor, final_factor);
}

SAMRAI_INLINE_KEYWORD void
QuatSysSolver::setLogging(bool logging)
{
  d_fac_ops->enableLogging(logging);
}

