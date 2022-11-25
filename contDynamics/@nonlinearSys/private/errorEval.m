function [perfInd] = errorEval(IHerrorActual,IHerrorAssume)
% errorEval - computes performance index of the linearization error
%
% Syntax:  
%    [perfInd] = errorEval(IHerrorActual,IHerrorAssume)
%
% Inputs:
%    IHerrorActual - cell array of actual linearizatuion errors
%    IHerrorAssume - assumed linearization error 
%
% Outputs:
%    perfInd - performance index
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Matthias Althoff
% Written: 15-January-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------


%compute performance index 
perfInd=max(edgeLength(IHerrorActual)./edgeLength(IHerrorAssume));


%------------- END OF CODE --------------