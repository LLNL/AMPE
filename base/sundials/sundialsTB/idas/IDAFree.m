function [] = IDAFree()
%IDAFree deallocates memory for the IDAS solver.
%
%   Usage:  IDAFree

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/08/21 17:38:42 $

mode = 40;
idm(mode);
