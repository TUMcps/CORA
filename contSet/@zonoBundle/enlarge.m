function zB = enlarge(zB,factor)
% enlarge - enlarges the generators of a zonotope bundle by a vector of
%    scaling factors
%
% Syntax:
%    zB = enlarge(zB,factor)
%
% Inputs:
%    zB - zonoBundle object
%    factor - vector of factors for the enlargement of each dimension
%
% Outputs:
%    zB - enlarged zonotope bundle
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    factor = [1.5;2];
%
%    res = enlarge(zB,factor);
%
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       20-November-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i=1:zB.parallelSets
    zB.Z{i} = enlarge(zB.Z{i},factor);
end

% ------------------------------ END OF CODE ------------------------------
