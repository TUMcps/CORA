function [Rfirst] = coreReach(obj,Rinit)
% coreReach - computes the reachable continuous set for the next time step
% without considering uncertain inputs
%
% Syntax:  
%    [Rfirst] = coreReach(obj,Rinit,options)
%
% Inputs:
%    Rinit - initial reachable set
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
% Written:      03-May-2011 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object structure
eAt=obj.taylor.eAt;
Rtrans=obj.taylor.Rtrans;

%first time step homogeneous solution
Rhom_tp=eAt*Rinit + Rtrans;

%write results to reachable set struct Rfirst
Rfirst=Rhom_tp;


%------------- END OF CODE --------------