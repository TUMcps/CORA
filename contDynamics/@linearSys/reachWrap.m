function [Rnext,options] = reachWrap(obj,R,options)
% reachWrap - classical reachability with wrapping effect
%
% Syntax:  
%    [Rnext,options] = reachWrap(obj,R,options)
%
% Inputs:
%    obj - linearSys object
%    R - reachable set of the previous time step
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      07-November-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%overall solution
R.ti=obj.taylor.eAt*R.ti+obj.taylor.RV+obj.taylor.Rtrans; 

%reduce zonotope
Rred=reduce(R.ti,'girard',options.zonotopeOrder);

%write results to reachable set struct Rnext
Rnext.tp=[];
Rnext.ti=Rred;



%------------- END OF CODE --------------