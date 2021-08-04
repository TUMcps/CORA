function [OGain, tComp] = observe_gain_PRadA(obj,options)
% observe_gain_PRadA - computes the gain for the guaranted state estimation
% approach according to Sec. 4.1 of [1].
%
% Syntax:  
%    [R,Rout] = observe_PRadA(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
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
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       12-Sep-2020
% Last update:   02-Jan-2021
%                25-Feb-2021
% Last revision: ---


%------------- BEGIN CODE --------------

tic

% store full output matrix and noise set
tmpC = obj.C;
V = options.V;

% obtain system dimension and nr of outputs
nrOfOutputs = size(tmpC,1);

% find optimal gain for each output
for i = 1:nrOfOutputs
    
    % extract i-th output
    obj.C = tmpC(i,:);
    options.V = box(project(V,i));
    
    % PRadA computes for all inputs simultaneously; using only a single
    % output is a special case
    % The approach in Sec. 4.1 of [1] can be seen as a special case of [2] 
    OGain(:,i) = observe_gain_PRadB(obj,options);
end

% restore old C matrix since system object makes only shallow copies
obj.C = tmpC;

% computation time
tComp = toc;

%------------- END OF CODE --------------