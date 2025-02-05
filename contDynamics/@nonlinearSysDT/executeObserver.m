function [R,tcomp] = executeObserver(nlnsysDT,params,options)
% executeObserver - calls the appropriate observer
%
% Syntax:
%    [R,tcomp] = executeObserver(nlnsysDT,params,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    params - model parameters for the guaranteed state estimation
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
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
switch options.alg
    case 'VolMin-A'
        [R,tcomp] = priv_observe_volMinA(nlnsysDT,params,options);
    case 'VolMin-B'
        [R,tcomp] = priv_observe_volMinB(nlnsysDT,params,options);
    case 'FRad-A'
        [R,tcomp] = priv_observe_FRadA(nlnsysDT,params,options);
    case 'FRad-B'
        [R,tcomp] = priv_observe_FRadB(nlnsysDT,params,options);
    case 'FRad-C'
        [R,tcomp] = priv_observe_FRadC(nlnsysDT,params,options);
end

% ------------------------------ END OF CODE ------------------------------
