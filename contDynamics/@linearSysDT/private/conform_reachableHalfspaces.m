function [N, d, E] = conform_reachableHalfspaces(obj,params,options)
% conform_reachableHalfspaces - computes reachable halfspaces according to 
% [1]. This is not realized by first computing reach() becasue only the 
% difference to the nominal solution is required.
%
% Syntax:
%    res = conform_reachableHalfspaces(obj,params,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    N - normal vectors of reachable halfsaces
%    d - offset of reachable halfspaces
%    E - precomputed matrices (E_i = C*A^i; see eq. (17) in [1])
%
%    [1] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       29-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load variables

A = obj.A;
C = obj.C;

% maximum number of timeSteps
maxNrOfTimeSteps = ceil(params.tFinal/obj.dt); 

%% Compute E_i, sum_E and reachable sets, see (eq. (15), (17) in [1])
% indices are shifted by one compared to [1]
% initialize cell array of E values
E = cell(maxNrOfTimeSteps,1);
N = cell(maxNrOfTimeSteps,1);
d = cell(maxNrOfTimeSteps,1);
E{1} = C;
Adi = A;

% deviation from center of initial states
X0 = params.R0conf;

% initialize reachable sets
R = obj.C*X0 + params.V; 
Radd = 0;

% Compute halfspace normals
Rtmp = halfspace(R); 
N{1} = Rtmp.halfspace.H;
d{1} = Rtmp.halfspace.K;

% loop over remaining halfspaces
for k = 2:maxNrOfTimeSteps+1
    E{k} = C*Adi; % E_i
    Radd = Radd + E{k-1}*params.W; % newly added partial reachable set
    R = E{k}*X0 + Radd + params.V; % reachable set
    Adi = A*Adi; % update i-th power of A
    % Compute halfspace normals
    R = halfspace(R); 
    N{k} = R.halfspace.H;
    d{k} = R.halfspace.K;
end


% ------------------------------ END OF CODE ------------------------------
