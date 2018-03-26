function status = CVodeAdjInit(steps, interp)
%CVodeAdjInit allocates and initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjInit(STEPS, INTEPR) 
%
%   STEPS    specifies the (maximum) number of integration steps between two 
%            consecutive check points.
%   INTERP   Specifies the type of interpolation used for estimating the forward 
%            solution during the backward integration phase. INTERP should be
%            'Hermite', indicating cubic Hermite interpolation, or 'Polynomial',
%            indicating variable order polynomial interpolation.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/12/05 21:58:17 $

mode = 4;

if nargin ~= 2
  error('Wrong number of input arguments');
end

status = cvm(mode,steps,interp);
