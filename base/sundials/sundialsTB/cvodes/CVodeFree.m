function CVodeFree()
%CVodeFree deallocates memory for the CVODES solver.
%
%   Usage:  CVodeFree

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.4 $Date: 2007/05/11 18:51:31 $

mode = 40;
cvm(mode);
