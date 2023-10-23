function cZ = minkDiff(cZ1,S,varargin)
% minkDiff - compute the Minkowski difference of two constrained zonotopes:
%         cZ1 - cZ2 = cZ <-> cZ + cZ2 \subseteq cZ1
%
% Syntax:
%    cZ = minkDiff(cZ1,S)
%    cZ = minkDiff(cZ1,S,type)
%
% Inputs:
%    cZ1 - conZonotope object
%    S - conZonotope object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'inner')
%
% Outputs:
%    cZ - conZonotope object after Minkowski difference
%
% Example: 
%    Z = [0 4 0 0;0 0 4 0];
%    A = [1 1 1.5]; b = -0.5;
%    cZ1 = conZonotope(Z,A,b);
%
%    Z = [1 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ2 = conZonotope(Z,A,b);
%
%    cZ = minkDiff(cZ1,cZ2);
%
%    figure; hold on;
%    plot(cZ1);
%    plot(cZ2,[1,2],'r');
%    plot(cZ,[1,2],'g');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   09-November-2022 (MW, rename 'minkDiff')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different algorithms for different set representations
if isnumeric(S)
   cZ = cZ1 + (-S); 
   
elseif isa(S,'zonotope') || isa(S,'interval')
    
    % convert to zonotope
    Z = zonotope(S);
    
    % compute Minkowski difference according to Theorem 1 in [1]
    c = center(Z);
    G = generators(Z);

    cZ = cZ1 + (-c);

    for i = 1:size(G,2)
        cZ = and_(cZ + G(:,i),cZ + (-G(:,i)),'exact');
    end
    
elseif isa(S,'conZonotope') || isa(S,'polytope') || ...
       isa(S,'zonoBundle')
   
    % parse input arguments
    type = 'exact';
    if nargin > 2 && ~isempty(varargin{1})
        type = varargin{1};
    end
    
    % compute exact result or fast inner-approximation
    if strcmp(type,'inner')
       cZ = minkDiff(S,zonotope(S));
       return;
    end
    
    % compute exact result according to Lemma 1 in [1]
    V = vertices(S);  cZ = cZ1 - V(:,1);
    
    for i = 2:size(V,2)
       cZ = and_(cZ,cZ1 - V(:,i),'exact'); 
    end
    
else
    
    % parse input arguments
    type = 'exact';
    if nargin > 2 && ~isempty(varargin{1})
        type = varargin{1};
    end
    
    if strcmp(type,'exact')
        throw(CORAerror('CORA:noExactAlg',cZ1,S));
    end
    
    % compute inner-approximation by enclosing second set with zonotope
    cZ = minkDiff(cZ1,zonotope(S));
end

% ------------------------------ END OF CODE ------------------------------
