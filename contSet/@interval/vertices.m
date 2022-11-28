function V = vertices(I)
% vertices - Computes vertices of an interval object
%
% Syntax:  
%    V = vertices(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    I = interval([1; -1], [2; 1]);
%    V = vertices(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/vertices

% Author:       Matthias Althoff
% Written:      24-July-2006 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_vertices('interval',I);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    V = vars{1}; return
end

% compute matrix with all possible generator combinations
one_pm = [1 -1];
for i=1:length(I.inf)-1
    one_pm = [ones(1,2^i) -ones(1,2^i); one_pm one_pm];
end

Iextended = [ones(1,size(one_pm,2));one_pm];

% convert to zonotope 
zono = zonotope(I);
Z = zono.Z;

% obtain vertices
V = zeros(size(Z,1),size(one_pm,2));
for i = 1:size(one_pm,2)
    V(:,i) = Z*Iextended(:,i);
end

%------------- END OF CODE --------------