function Rtp = linReach(obj,Rinit,options)
% linReach - computes the reachable set after linearization
%
% Syntax:  
%    [Rtp] = linReach(obj,Rinit,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    Rinit - initial reachable set
%    options - options struct
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
if options.tensorOrder == 3
    Verror = linError_thirdOrder(obj, options, Rdelta); 
elseif options.tensorOrder == 2
    Verror = linError_mixed_noInt(obj, options, Rdelta);   
else
    error("Tensor Order not implemented!");
end


%add interval of actual error
Rtp=Rtp+Verror;


%------------- END OF CODE --------------