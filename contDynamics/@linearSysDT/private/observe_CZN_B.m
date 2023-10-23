function [R,tcomp] = observe_CZN_B(obj,options)
% observe_CZN_B - computes the guaranteed state estimation approach
% from [1], [2], where the intersection method is from [1] and the
% propagation is described in [2].
%
%
% Syntax:
%    [R,tcomp] = observe_CZN_B(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
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
% Example: 
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
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% width of strips
sigma = supremum(abs(interval(options.V)));

%% Compute gain for intersection
F = generators(options.V);
E = generators(options.W);
Qv = F*F';
Qw = E*E';
% generators
G = options.R0.G;  
Pbar = G*G';
Rbar = obj.A*Pbar*obj.A' +Qw; 
S = obj.C*Rbar*obj.C' + Qv;
L = Rbar*obj.C'; 
Lambda = L / S; %L*inv(S); 

% Intersection of strips according to Theorem 6.3 of [1]
y = options.y(:,1);
Rnext.tp = intersectStrip(options.R0,obj.C,sigma,y,Lambda);

% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    
    % Prediction, part of eq. (33) of [2]
    Rnext.tp = obj.A*Rnext.tp + obj.B*options.uTransVec(:,k) + options.W;
    
    %% Update gain for intersection
    % generators
    G = Rnext.tp.G;  
    Pbar = G*G';
    Rbar = obj.A*Pbar*obj.A' +Qw; 
    S = obj.C*Rbar*obj.C' + Qv;
    L = Rbar*obj.C'; 
    Lambda = L / S; % L*inv(S); 
    
    % Intersection according to Theorem 6.3 of [1]
    y = options.y(:,k+1);
    Rnext.tp = intersectStrip(Rnext.tp,obj.C,sigma,y,Lambda);
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp;
end
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
