function [OGain, tComp] = observe_gain_PRadA(linsysDT,params,options)
% observe_gain_PRadA - computes the gain for the guaranteed state estimation
% approach according to Sec. 4.1 of [1].
%
% Syntax:
%    [R,Rout] = observe_gain_PRadA(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    OGain - observer gain
%    tComp - computation time
%
% Reference:
%    [1] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotopic guaranteed state estimation for
%        uncertain systems. Automatica, 49(11):3418–3424, 2013.
%    [2] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotope-based set-membership estimation for
%        multi-output uncertain systems. In Proc. of the IEEE
%        International Symposium on Intelligent Control (ISIC),
%        pages 212–217, 2013.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2020
% Last update:   02-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic;

% store full output matrix and noise set
tmpC = linsysDT.C;
V = params.V;

% obtain system dimension and nr of outputs
nrOfOutputs = linsysDT.nrOfOutputs;

% find optimal gain for each output
for i = 1:nrOfOutputs
    
    % extract i-th output
    linsysDT.C = tmpC(i,:);
    params.V = box(project(V,i));
    
    % PRadA computes for all inputs simultaneously; using only a single
    % output is a special case
    % The approach in Sec. 4.1 of [1] can be seen as a special case of [2] 
    OGain(:,i) = observe_gain_PRadB(linsysDT,params,options);
end

% restore old C matrix since system object makes only shallow copies
linsysDT.C = tmpC;

% computation time
tComp = toc;

% ------------------------------ END OF CODE ------------------------------
