function [Rnext,options] = reachAlternative(obj,R,options)
% reachAlternative - computes the reachable continuous probabilistic set for one time 
% step of a linear interval system
%
% Syntax:  
%    [Rnext,options] = reachAlternative(obj,R,options)
%
% Inputs:
%    obj - linProbSys object
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
% Written:      27-September-2007 
% Last update:  09-September-2009
% Last revision:---

%------------- BEGIN CODE --------------

%next set (no input)
R.ti=obj.taylor.eAt*R.ti;

%bloating due to input
R.ti=R.ti+obj.taylor.Rinput+obj.taylor.pRinput;

%reduce
Rnext.ti=reduce(R.ti,'girard',options.zonotopeOrder);
Rnext.tp=[];

%------------- END OF CODE --------------