function [qZ1] = exactPlus(qZ1,qZ2)
% exactPlus - Adds two quadZonotopes by adding all generators of common
% generator factors, except for Grest. It has to be ensured from outside
% that the generator factors match
%
% Syntax:  
%    [qZ] = exactPlus(summand1,summand2)
%
% Inputs:
%    qZ1 - quadZonotope object
%    qZ2 - quadZonotope object
%
% Outputs:
%    qZ1 - quadZonotpe after Minkowsi addition
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      05-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%addition
qZ1.c = qZ1.c + qZ2.c;
qZ1.G = qZ1.G + qZ2.G;
qZ1.Gsquare = qZ1.Gsquare + qZ2.Gsquare;
qZ1.Gquad = qZ1.Gquad + qZ2.Gquad;

%concatenation
qZ1.Grest(:,end+1:end+length(qZ2.Grest(1,:))) = qZ2.Grest; 


%concatenation

%------------- END OF CODE --------------