function [Rnext,options] = initReach_adaptive(obj,options)
% initReach_adaptive - computes the reachable continuous set
%   for the first time step
%
% Syntax:
%    [Rnext,options] = initReach_adaptive(obj,options)
%
% Inputs:
%    obj - nonlinearSys object
%    options - struct containing the algorithm settings
%
% Outputs:
%    Rnext - reachable set
%    options - struct containing the algorithm settings
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       14-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[Rti,Rtp,~,options] = linReach_adaptive(obj,options,options.R);

% store the results
Rnext.tp = Rtp;
Rnext.ti = Rti;
Rnext.R0 = options.R;

% ------------------------------ END OF CODE ------------------------------
