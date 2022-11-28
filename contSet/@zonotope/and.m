function Z = and(Z,S,varargin)
% and - overloads & operator, computes the intersection of two zonotopes
%
% Syntax:  
%    Z = and(Z,S)
%    Z = and(Z,S,method)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%    method - (optional) algorithm used to compute the intersection
%               - 'conZonotope' (default)
%               - 'averaging'
%
% Outputs:
%    Z - zonotope object enclosing the intersection 
%
% Example: 
%    Z1 = zonotope([4 2 2;1 2 0]);
%    Z2 = zonotope([3 1 -1 1;3 1 2 0]);
%
%    res = Z1 & Z2;
%
%    figure; hold on;
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'b');
%    plot(res,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/and

% Author:        Matthias Althoff, Niklas Kochdumper, Amr Alanwar
% Written:       29-June-2009
% Last update:   02-Sep-2019 (rename intersection -> and)
%                11-Nov-2019 (NK: added algorithm for general case)
%                12-Feb-2020 (Amr: adding averaging option)
% Last revision: ---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_and('zonotope',Z,S,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    Z = vars{1}; return
else
    Z = vars{1}; S = vars{2}; method = vars{3};
end


% quick check: simpler function for intervals
if isInterval(Z) && isInterval(S)
    % conversion to intervals exact
    Z = zonotope(interval(Z) & interval(S)); return
end

% special algorithm for two parallelotopes with the same center
if isa(S,'levelSet') || isa(S,'conZonotope') || ...
   isa(S,'zonoBundle') || isa(S,'conPolyZono')
    
    Z = S & Z;
    
elseif strcmp(method,'conZonotope')
    
    % convert sets to constrained zonotopes
    Z = conZonotope(Z);
    
    if ~isa(S,'halfspace') && ~isa(S,'conHyperplane')
        S = conZonotope(S);
    end
    
    % compute intersection
    Z = Z & S;
    
    % conclose resulting constrained zonotope by a zonotope
    Z = zonotope(Z);

elseif strcmp(method,'averaging')
    
    list{1} = Z;
    list{2} = S;
    Z = andAveraging(list);
    
else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',Z,S));
end

%------------- END OF CODE --------------