function [] = my_install_STB
%
% INSTALL_STB Interactive compilation and installtion of sundialsTB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.15 $Date: 2007/02/05 20:20:33 $

% MEX compiler command
% --------------------

mexcompiler = 'mex -g';

% Location of sundialsTB and top of sundials source tree
% ------------------------------------------------------

stb = pwd;
cd('..');
sun = pwd;
cd(stb);

% Should we enable parallel support?
% ----------------------------------

par = true;
if isempty(getenv('LAMHOME'))
    par = false;
end
if isempty(getenv('MPITB_ROOT'))
    par = false;
end
q = fullfile(sun,'src','nvec_par');
if ~exist(q, 'dir')
    par = false;
end

% Create sundials_config.h
% ------------------------

mkdir('sundials');
fi = fopen(fullfile('sundials','sundials_config.h'),'w');
fprintf(fi,'#define SUNDIALS_PACKAGE_VERSION "2.3.0"\n');
fprintf(fi,'#define SUNDIALS_DOUBLE_PRECISION 1\n');
fprintf(fi,'#define SUNDIALS_USE_GENERIC_MATH 1\n');
fprintf(fi,'#define SUNDIALS_EXPORT\n');
fclose(fi);

% Compile MEX file
% ----------------


cvm_ok = true;
kim_ok = true;
idm_ok = true;
%compile_CVM(mexcompiler,stb,sun,par);
%compile_IDM(mexcompiler,stb,sun,par);
compile_KIM(mexcompiler,stb,sun,par);

% Remove sundials_config.h
% ------------------------

rmdir('sundials','s');

fprintf('\nMEX files were successfully created.\n');

%---------------------------------------------------------------------------------
% compilation of cvm MEX file
%---------------------------------------------------------------------------------

function [] = compile_CVM(mexcompiler,stb,sun,par)

cvm_sources = {
    fullfile(stb,'cvodes','cvm','src','cvm.c')
    fullfile(stb,'cvodes','cvm','src','cvmWrap.c')
    fullfile(stb,'cvodes','cvm','src','cvmOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(cvm_sources)
    sources = sprintf('%s "%s"',sources,cvm_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

cvm_incdir = fullfile(stb,'cvodes','cvm','src'); % for cvm.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,cvm_incdir,nvm_incdir);

libraries = '';

% Add CVODES sources and header files

cvs_sources = {
    fullfile(sun,'src','cvodes','cvodes_band.c')
    fullfile(sun,'src','cvodes','cvodes_bandpre.c')
    fullfile(sun,'src','cvodes','cvodes_bbdpre.c')
    fullfile(sun,'src','cvodes','cvodes_direct.c')
    fullfile(sun,'src','cvodes','cvodes_dense.c')
    fullfile(sun,'src','cvodes','cvodes_diag.c')
    fullfile(sun,'src','cvodes','cvodea.c')
    fullfile(sun,'src','cvodes','cvodes.c')
    fullfile(sun,'src','cvodes','cvodes_io.c')
    fullfile(sun,'src','cvodes','cvodea_io.c')
    fullfile(sun,'src','cvodes','cvodes_spils.c')
    fullfile(sun,'src','cvodes','cvodes_spbcgs.c')
    fullfile(sun,'src','cvodes','cvodes_spgmr.c')
    fullfile(sun,'src','cvodes','cvodes_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };
for i=1:length(cvs_sources)
    sources = sprintf('%s "%s"',sources,cvs_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');      % for SUNDIALS exported headers
cvs_srcdir = fullfile(sun,'src','cvodes'); % for cvodes_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,cvs_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

cvm_dir = fullfile(stb,'cvodes','cvm');
cd(cvm_dir)
mex_cmd = sprintf('%s -O -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)

%---------------------------------------------------------------------------------
% compilation of idm MEX file
%---------------------------------------------------------------------------------

function [] = compile_IDM(mexcompiler,stb,sun,par)

idm_sources = {
    fullfile(stb,'idas','idm','src','idm.c')
    fullfile(stb,'idas','idm','src','idmWrap.c')
    fullfile(stb,'idas','idm','src','idmOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(idm_sources)
    sources = sprintf('%s "%s"',sources,idm_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

idm_incdir = fullfile(stb,'idas','idm','src');   % for idm.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,idm_incdir,nvm_incdir);

libraries = '';

% Add IDAS sources and header files

ids_sources = {
    fullfile(sun,'src','idas','idas_band.c')
    fullfile(sun,'src','idas','idas_bbdpre.c')
    fullfile(sun,'src','idas','idas_dense.c')
    fullfile(sun,'src','idas','idas_direct.c')
    fullfile(sun,'src','idas','idaa.c')
    fullfile(sun,'src','idas','idas.c')
    fullfile(sun,'src','idas','idas_ic.c')
    fullfile(sun,'src','idas','idas_io.c')
    fullfile(sun,'src','idas','idaa_io.c')
    fullfile(sun,'src','idas','idas_spils.c')
    fullfile(sun,'src','idas','idas_spbcgs.c')
    fullfile(sun,'src','idas','idas_spgmr.c')
    fullfile(sun,'src','idas','idas_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };
for i=1:length(ids_sources)
    sources = sprintf('%s "%s"',sources,ids_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');     % for SUNDIALS exported headers
ids_srcdir = fullfile(sun,'src','idas');  % for idas_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,ids_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

idm_dir = fullfile(stb,'idas','idm');
cd(idm_dir)
mex_cmd = sprintf('%s -O -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)


%---------------------------------------------------------------------------------
% compilation of KINSOL MEX file
%---------------------------------------------------------------------------------

function [] = compile_KIM(mexcompiler,stb,sun,par)

kim_sources = {
    fullfile(stb,'kinsol','kim','src','kim.c')
    fullfile(stb,'kinsol','kim','src','kimWrap.c')
    fullfile(stb,'kinsol','kim','src','kimOpts.c')
              };
if par
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_parallel.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
else
    nvm_sources = {
        fullfile(stb,'nvector','src','nvm_serial.c')
        fullfile(stb,'nvector','src','nvm_ops.c')
                  };
end
sources = '';
for i=1:length(kim_sources)
    sources = sprintf('%s "%s"',sources,kim_sources{i});
end
for i=1:length(nvm_sources)
    sources = sprintf('%s "%s"',sources,nvm_sources{i});
end

kim_incdir = fullfile(stb,'kinsol','kim','src'); % for kim.h
nvm_incdir = fullfile(stb,'nvector','src');      % for nvm.h
includes = sprintf('-I"%s" -I"%s" -I"%s"',stb,kim_incdir,nvm_incdir);

libraries = '';

% Add KINSOL sources and header files

kin_sources = {
    fullfile(sun,'src','kinsol','kinsol_band.c')
    fullfile(sun,'src','kinsol','kinsol_bbdpre.c')
    fullfile(sun,'src','kinsol','kinsol_dense.c')
    fullfile(sun,'src','kinsol','kinsol_direct.c')
    fullfile(sun,'src','kinsol','kinsol.c')
    fullfile(sun,'src','kinsol','kinsol_io.c')
    fullfile(sun,'src','kinsol','kinsol_spils.c')
    fullfile(sun,'src','kinsol','kinsol_spbcgs.c')
    fullfile(sun,'src','kinsol','kinsol_spgmr.c')
    fullfile(sun,'src','kinsol','kinsol_sptfqmr.c')
              };
shr_sources = {
    fullfile(sun,'src','sundials','sundials_band.c')
    fullfile(sun,'src','sundials','sundials_dense.c')
    fullfile(sun,'src','sundials','sundials_iterative.c')
    fullfile(sun,'src','sundials','sundials_nvector.c')
    fullfile(sun,'src','sundials','sundials_direct.c')
    fullfile(sun,'src','sundials','sundials_spbcgs.c')
    fullfile(sun,'src','sundials','sundials_spgmr.c')
    fullfile(sun,'src','sundials','sundials_sptfqmr.c')
    fullfile(sun,'src','sundials','sundials_math.c')
              };

for i=1:length(kin_sources)
    sources = sprintf('%s "%s"',sources,kin_sources{i});
end
for i=1:length(shr_sources)
    sources = sprintf('%s "%s"',sources,shr_sources{i});
end

sun_incdir = fullfile(sun,'include');      % for SUNDIALS exported headers
kin_srcdir = fullfile(sun,'src','kinsol'); % for kinsol_impl.h
includes = sprintf('%s -I"%s" -I"%s"',includes,sun_incdir,kin_srcdir);

% Add NVEC_SER sources and header files

nvs_sources = fullfile(sun,'src','nvec_ser','nvector_serial.c');
sources = sprintf('%s "%s"',sources, nvs_sources);

if par
  
  % Add NVEC_PAR sources and header files

  nvp_sources = fullfile(sun,'src','nvec_par','nvector_parallel.c');
  sources = sprintf('%s "%s"',sources, nvp_sources);

% Add LAM headers and libraries

  lam = getenv('LAMHOME');
  lam_incdir = fullfile(lam, 'include');
  lam_libdir = fullfile(lam, 'lib');
  includes = sprintf('%s -I"%s"',includes,lam_incdir);
  libraries = sprintf('%s -L"%s" -lmpi -llam -lutil',libraries,lam_libdir);

end

% Create MEX file

kim_dir = fullfile(stb, 'kinsol', 'kim');
cd(kim_dir)
mex_cmd = sprintf('%s -v %s %s %s', mexcompiler, includes, sources, libraries);
disp(mex_cmd);
eval(mex_cmd);

% Move back to sundialsTB

cd(stb)

