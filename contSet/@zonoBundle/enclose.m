function zB = enclose(varargin)
% enclose - Generates a zonoBundle that encloses a zonoBundle and its
%           linear transformation (see Proposition 5 in [1])
%
% Syntax:  
%    zB = enclose(zB1,zB2)
%    zB = enclose(zB1,M,zBplus)
%
% Inputs:
%    zB1 - first zonoBundle object
%    zB2 - second zonoBundle object, satisfying zB2 = (M * zB1) + zBplus
%    M - matrix for the linear transformation
%    zBplus - zonoBundle object added to the linear transformation
%
% Outputs:
%    zB - zonoBundle that encloses zB1 and zB2
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Matthias Althoff
% Written:      10-November-2010 
% Last update:  25-January-2016
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 2
    Zbundle1 = varargin{1};
    Zbundle2 = varargin{2};
else
    Zbundle1 = varargin{1};
    M = varargin{2};
    zBplus = varargin{3};
    
    Zbundle2 = (M*Zbundle1) + zBplus;
end

% compute set
zB = Zbundle1;

if isa(Zbundle2,'zonoBundle')
    %compute enclosure for each zonotope pair
    for i=1:Zbundle1.parallelSets
        zB.Z{i}=enclose(Zbundle1.Z{i},Zbundle2.Z{i});
    end
elseif isa(Zbundle2,'zonotope')
    %compute enclosure for each zonotope
    for i=1:Zbundle1.parallelSets
        zB.Z{i}=enclose(Zbundle1.Z{i},Zbundle2);
    end
end


%------------- END OF CODE --------------