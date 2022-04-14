function [R,tcomp] = executeObserver(obj,options)
% executeObserver - calls the appropriate observer
%
% Syntax:  
%    [R,tcomp] = executeObserver(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
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

% Author:       Matthias Althoff
% Written:      20-Mar-2020
% Last update:  25-Feb-2021
%               14-Jun-2021
% Last revision:---


%------------- BEGIN CODE --------------

% decide which observer to execute by options.alg
if strcmp(options.alg,'VolMin-A') 
    [R,tcomp] = observe_volMinA(obj, options);
elseif strcmp(options.alg,'VolMin-B') 
    [R,tcomp] = observe_volMinB(obj, options);
elseif strcmp(options.alg,'FRad-A') 
    [R,tcomp] = observe_FRadA(obj, options);
elseif strcmp(options.alg,'FRad-B') 
    [R,tcomp] = observe_FRadB(obj, options);
elseif strcmp(options.alg,'PRad-A') 
    [R,tcomp] = observe_PRadA(obj, options);
elseif strcmp(options.alg,'PRad-B') 
    [R,tcomp] = observe_PRadB(obj, options);
elseif strcmp(options.alg,'PRad-C') 
    [R,tcomp] = observe_PRadC(obj, options);
elseif strcmp(options.alg,'FRad-C') 
    [R,tcomp] = observe_FRadC(obj, options);
elseif strcmp(options.alg,'PRad-D') 
    [R,tcomp] = observe_PRadD(obj, options);
elseif strcmp(options.alg,'PRad-E') 
    [R,tcomp] = observe_PRadE(obj, options);
elseif strcmp(options.alg,'Nom-G') 
    [R,tcomp] = observe_NomG(obj, options);
elseif strcmp(options.alg,'Hinf-G') 
    [R,tcomp] = observe_HinfG(obj, options);
elseif strcmp(options.alg,'ESO-A')
    [R,tcomp] = observe_ESO_A(obj, options);
elseif strcmp(options.alg,'ESO-B')
    [R,tcomp] = observe_ESO_B(obj, options);
elseif strcmp(options.alg,'ESO-C')
    [R,tcomp] = observe_ESO_C(obj, options);    
elseif strcmp(options.alg,'ESO-D')
    [R,tcomp] = observe_ESO_D(obj, options);
elseif strcmp(options.alg,'CZN-A')
    [R,tcomp] = observe_CZN_A(obj, options);  
elseif strcmp(options.alg,'CZN-B')
    [R,tcomp] = observe_CZN_B(obj, options);  
elseif strcmp(options.alg,'ROPO')
    [R,tcomp] = observe_ROPO(obj,options);
end

%------------- END OF CODE --------------

