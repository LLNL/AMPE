function status = IDAAdjReInit()
%IDAAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: IDAAdjReInit
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:18 $

mode = 14;

status = idm(mode);
