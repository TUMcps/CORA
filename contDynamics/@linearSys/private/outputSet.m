function [Rout,Rout_tp] = outputSet(C,D,k,R,options)
% outputSet - calculates output set based on output equation given by
%  y = Cx + Du + k and sets for x (R) and u(options.U|uTrans)
% including removal of empty generators
%
% Syntax:  
%    [Rout,Rout_tp] = outputSet(C,D,k,R,options)
%
% Inputs:
%    (matrices from y = Cx + Du + k)
%    R       - reachable set of current step
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rout    - output set of time intervals
%    Rout_tp - output set of time points
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      12-Aug-2019
% Last update:  20-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

Z    = R.ti; % time interval solution
Z_tp = R.tp; % time point solution


if isfield(options,'saveOrder')
    Rout = reduce(zonotope(C*Z) + D * (options.uTrans + options.U) + k,...
        options.reductionTechnique,options.saveOrder);
    Rout_tp = reduce(zonotope(C*Z_tp) + D * (options.uTrans + options.U) + k,...
        options.reductionTechnique,options.saveOrder);
else
    if any(D)
        Rout = C*Z + D * (options.uTrans + options.U) + k;
        Rout_tp = C*Z_tp + D * (options.uTrans + options.U) + k;
    elseif isscalar(C) && C == 1 && ~any(k)
        % speed up for systems without output equation
        Rout = Z;
        Rout_tp = Z_tp;
    else
        Rout = C*Z + k;
        Rout_tp = C*Z_tp + k;
    end
end


end



%------------- END OF CODE --------------

