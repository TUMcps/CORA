function zB = shrink(zB,filterLength)
% shrink - shrinks each zonotope in a zonotope bundle; however: in total,
%    the size of the zonotope bundle will increase
%
% Syntax:  
%    zB = shrink(zB,filterLength)
%
% Inputs:
%    zB - zonoBundle object
%    filterLength - ???
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    ---
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      01-December-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% compute over-approximative parallelotope for each zonotope
Zred = cell(zB.parallelSets,1);
for i=1:zB.parallelSets
    Zred{i} = reduce(zB.Z{i},'methC',1,filterLength);
end

% intersect all parallelotopes
P = mptPolytope(Zred{1});
for i=2:zB.parallelSets
    P = and_(P,mptPolytope(Zred{i}),'exact');
end

% check if polytope is empty
if ~isempty(P)
    % over-approximate polytope by zonotopes
    for i=1:zB.parallelSets
        Zmat = Zred{i}.Z;
        zB.Z{i} = parallelotope(P,Zmat(:,2:end));
    end
else
    zB = [];
end

%------------- END OF CODE --------------