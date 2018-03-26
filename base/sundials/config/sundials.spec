Summary: suite of differential and nonlinear equation solvers
Name: sundials
Version: 2.3.0
Release: current
License: BSD-like
Group: Development/Libraries
Source: http://www.llnl.gov/casc/sundials/download/%{name}-%{version}.tgz
URL: http://www.llnl.gov/casc/sundials/
BuildRoot: %_tmppath/%name-%version-%release-root-%(%__id_u -n)

%description
SUNDIALS is a SUite of Non-linear DIfferential/ALgebraic equation Solvers
for use in writing mathematical software.

SUNDIALS was implemented with the goal of providing robust time integrators and 
nonlinear solvers that can easily be incorporated into existing simulation codes. 
SUNDIALS consists of the folowing solvers:
CVODE 	solves initial value problems for ordinary differential equation (ODE) systems.
CVODES 	solves ODE systems and includes sensitivity analysis capabilities (forward and adjoint).
IDA 	solves initial value problems for differential-algebraic equation (DAE) systems.
KINSOL 	solves nonlinear algebraic systems.

%prep
%setup -q 

%build
%configure --enable-shared --enable-static 
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc INSTALL_NOTES LICENSE README
%_libdir
%_includedir
%_bindir

%changelog
* Fri Jul 28 2006 Radu Serban <radu@llnl.gov>
- Changed summary and description
- Added bindir (for sundials-config)
- Changed tarball extension to 'tgz'
* Thu Jul 27 2006 John Pye <john.pye@student.unsw.edu.au> 2.3.0-0
- First RPM spec created

