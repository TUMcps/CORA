function [Rfirst,options] = initReach(obj,Rinit,options)
% reach - computes the reachable continuous set for the first time step
%
% Syntax:  
%    [Rnext] = reach(obj,R,options)
%
% Inputs:
%    obj - zeroDynSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%    Rfirst - first reachable set 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      02-November-2007 
% Last update:  06-July-2009
%               18-March-2017
% Last revision: ---

%------------- BEGIN CODE --------------

Rfirst.tp=Rinit;
Rfirst.ti=Rinit;

%------------- END OF CODE --------------