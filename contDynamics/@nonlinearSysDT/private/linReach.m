function Rtp = linReach(obj,Rinit,options)
% linReach - computes the reachable set after linearazation
%
% Syntax:  
%    [Rtp] = linReach(obj,Rinit,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    options - options struct
%    Rinit - initial reachable set
%
% Outputs:
%    Rtp - resulting reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

% linearize nonlinear system
[obj,A_lin,U] = linearize(obj,Rinit,options); 

%translate Rinit by linearization point
Rdelta = Rinit + (-obj.linError.p.x);

% compute reachable set of linearized system
Rtp = A_lin*Rdelta + U;

% obtain linearization error
if options.tensorOrder > 2
    Verror = linError_thirdOrder(obj, options, Rdelta); 
else
    Verror = linError_mixed_noInt(obj, options, Rdelta);   
end


%add interval of actual error
Rtp=Rtp+Verror;


%------------- END OF CODE --------------