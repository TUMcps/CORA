function [R,tcomp] = priv_observe_CZN_B(linsysDT,params,options)
% priv_observe_CZN_B - computes the guaranteed state estimation approach
%    from [1], [2], where the intersection method is from [1] and the
%    propagation is described in [2].
%
% Syntax:
%    [R,tcomp] = priv_observe_CZN_B(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] Amr Alanwar, Victor Gassmann, Xingkang He, Hazem Said,
%        Henrik Sandberg, Karl Henrik Johansson, and Matthias
%        Althoff. Privacy preserving set-based estimation using
%        partially homomorphic encryption. arXiV.org.
%    [2] J. K. Scott, D. M. Raimondo, G. R. Marseglia, and R. D.
%        Braatz. Constrained zonotopes: A new tool for set-based
%        estimation and fault detection. Automatica, 69:126–136,
%        2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic;

% time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% width of strips
sigma = supremum(abs(interval(params.V)));

%% Compute gain for intersection
F = generators(params.V);
E = generators(params.W);
Qv = F*F';
Qw = E*E';
% generators
G = params.R0.G;  
Pbar = G*G';
Rbar = linsysDT.A*Pbar*linsysDT.A' +Qw; 
S = linsysDT.C*Rbar*linsysDT.C' + Qv;
L = Rbar*linsysDT.C'; 
Lambda = L / S; %L*inv(S); 

% Intersection of strips according to Theorem 6.3 of [1]
y = params.y(:,1);
Rnext.tp = intersectStrip(params.R0,linsysDT.C,sigma,y,Lambda);

% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    % Prediction, part of eq. (33) of [2]
    Rnext.tp = linsysDT.A*Rnext.tp + linsysDT.B*params.uTransVec(:,k) + params.W;
    
    %% Update gain for intersection
    % generators
    G = Rnext.tp.G;  
    Pbar = G*G';
    Rbar = linsysDT.A*Pbar*linsysDT.A' +Qw; 
    S = linsysDT.C*Rbar*linsysDT.C' + Qv;
    L = Rbar*linsysDT.C'; 
    Lambda = L / S; % L*inv(S); 
    
    % Intersection according to Theorem 6.3 of [1]
    y = params.y(:,k+1);
    Rnext.tp = intersectStrip(Rnext.tp,linsysDT.C,sigma,y,Lambda);
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp;
end
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
