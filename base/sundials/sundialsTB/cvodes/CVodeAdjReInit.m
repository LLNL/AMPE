function status = CVodeAdjReInit()
%CVodeAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjReInit
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:17 $

mode = 14;

status = cvm(mode);
