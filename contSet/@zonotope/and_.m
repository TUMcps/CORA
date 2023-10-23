function Z = and_(Z,S,method,varargin)
% and_ - overloads & operator, computes the intersection of two zonotopes
%
% Syntax:
%    Z = and_(Z,S)
%    Z = and_(Z,S,method)
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
% See also: contSet/and, conPolyZono/and_

% Authors:       Matthias Althoff, Niklas Kochdumper, Amr Alanwar
% Written:       29-June-2009
% Last update:   02-September-2019 (rename intersection -> and)
%                11-November-2019 (NK, added algorithm for general case)
%                12-February-2020 (AA, adding averaging option)
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

% quick check: simpler function for intervals
if representsa_(Z,'interval',eps) && representsa_(S,'interval',eps)
    % conversion to intervals exact
    Z = zonotope(and_(interval(Z),interval(S),'exact')); return
end

% special algorithm for two parallelotopes with the same center
if isa(S,'levelSet') || isa(S,'conZonotope') || ...
   isa(S,'zonoBundle') || isa(S,'conPolyZono')
    
    Z = and_(S,Z,'exact');
    
elseif strcmp(method,'conZonotope')
    
    % convert sets to constrained zonotopes
    Z = conZonotope(Z);
    
    if ~isa(S,'halfspace') && ~isa(S,'conHyperplane')
        S = conZonotope(S);
    end
    
    % compute intersection
    Z = and_(Z,S,'exact');
    
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

% ------------------------------ END OF CODE ------------------------------
