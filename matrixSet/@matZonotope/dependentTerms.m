function [zSq,zH] = dependentTerms(obj,r)
% dependentTerms - computes exact Taylor terms of a matrix zonotope square and
% a matrix zonotope exponential
%
% These different tasks are computed in one m-file to save computation
% time: the for loop has to be executed only once and help functions do not
% have to be called so often
%
% Syntax:  
%    [zSq,zH] = dependentTerms(obj,r)
%
% Inputs:
%    obj - linear interval system (linIntSys) object
%    r - time step increment
%
% Outputs:
%    zSq - exact square matrix
%    zH - exact Taylor terms up to second order
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      24-September-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object structure
C=obj.center;
G=obj.generator;
dim=obj.dim;
gens=obj.gens;


%square computation--------------------------------------------------------
%new center
sqC = C^2*r^2;
for i=1:gens
    sqC = sqC+(0.5*G{i}^2)*r^2;
end

%get generators
sqG = [];
%1st set of generators
for i=1:gens
    sqG{end+1}=(C*G{i} + G{i}*C)*r^2;
end
%2nd set of generators
for i=1:gens
    sqG{end+1}=0.5*G{i}^2*r^2;
end
%get indices for 3rd set of generators
if (gens>=2)
    ind = combinator(gens,2,'c');
    for i=1:length(ind(:,1))
        sqG{end+1}=(G{ind(i,1)}*G{ind(i,2)} + G{ind(i,2)}*G{ind(i,1)})*r^2;
    end
end
%--------------------------------------------------------------------------


%H computation-------------------------------------------------------------
%new center
HC = eye(dim) + C*r + 0.5*sqC;

%get generators
HG = [];
%1st set of generators
for i=1:gens
    HG{end+1}= G{i}*r + 0.5*sqG{i};
end
%remaining set of generators
for i=(gens+1):length(sqG)
    HG{i}=0.5*sqG{i};
end
%--------------------------------------------------------------------------

%write as matrix zonotopes
zSq=matZonotope(sqC,sqG);
zH=matZonotope(HC,HG);

%------------- END OF CODE --------------