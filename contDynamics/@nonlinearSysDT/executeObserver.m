function [R,tcomp] = executeObserver(obj,options)
% executeObserver - calls the appropriate observer
%
% Syntax:
%    [R,tcomp] = executeObserver(obj,options)
%
% Inputs:
%    obj - discrete-time nonlinear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% decide which observer to execute by options.alg
if strcmp(options.alg,'VolMin-A') 
    [R,tcomp] = observe_volMinA(obj, options);
elseif strcmp(options.alg,'VolMin-B') 
    [R,tcomp] = observe_volMinB(obj, options);
elseif strcmp(options.alg,'FRad-A') 
    [R,tcomp] = observe_FRadA(obj, options);
elseif strcmp(options.alg,'FRad-B') 
    [R,tcomp] = observe_FRadB(obj, options);
elseif strcmp(options.alg,'FRad-C') 
    [R,tcomp] = observe_FRadC(obj, options);
end

% ------------------------------ END OF CODE ------------------------------
