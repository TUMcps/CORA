function Z = or(Z1, Zcell)
% or - Computes union of zonotopes; works best when directions of
% generators are similar
%
% Syntax:  
%    Z = or(Z1, Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Zcell - cell array of zonotope objects
%
% Outputs:
%   Z - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      20-September-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%init
Zmat = [];

%dimension
dim = length(center(Z1));

%obtain minimum number of generators every zonotope has
minNrOfVecs = length(Z1.Z(1,:));
for iSet = 1:length(Zcell)
    minNrOfVecs = min(minNrOfVecs, length(Zcell{iSet}.Z(1,:)));
end

%obtain Zcut
Zcut{1} = Z1.Z(:, 1 : minNrOfVecs);
for iSet = 1:length(Zcell)
    Zcut{iSet+1} = Zcell{iSet}.Z(:, 1 : minNrOfVecs);
end

%obtain Zadd
Zadd{1} = Z1.Z(:, minNrOfVecs+1 : end);
for iSet = 1:length(Zcell)
    Zadd{iSet+1} = Zcell{iSet}.Z(:, minNrOfVecs+1 : end);
end

%compute vertex sets for each set of generators
for iGen = 1:minNrOfVecs
    v(:,1) = Zcut{1}(:,iGen);
    for iSet = 2:length(Zcut)
        v_new = Zcut{iSet}(:,iGen);
        %if direction is correct
        if v(:,1).'*v_new > 0
            v(:, iSet) = v_new;
        %flip direction
        else
            v(:, iSet) = - v_new;
        end
    end
    
    %generate vertex object
    V = vertices(v);
    
    %compute enclosing zonotope
    Z_encl = zonotope(V);
    
    %concatenate enclosing zonotopes
    Zmat(:,end+1 : end+dim+1) = Z_encl.Z;
end

%create enclosing zonotope
Z = zonotope(Zmat);



%------------- END OF CODE --------------