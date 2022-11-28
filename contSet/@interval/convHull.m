function res = convHull(I,varargin)
% convHull - computes an enclosure for the convex hull of two intervals
%
% Syntax:  
%    res = convHull(I,S)
%
% Inputs:
%    I - interval object
%    S - contSet object
%
% Outputs:
%    res - interval enclosing the convex hull
%
% Example:
%    I = interval([-1;-2],[3;4]);
%    Z = zonotope(ones(2,1),[1 0.5; -0.2 -2]);
%    res = convHull(I,Z);
%
%    figure; hold on;
%    plot(I,[1,2],'b');
%    plot(Z,[1,2],'b');
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:        Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
%                12-March-2021 (MW, add empty case)
% Last revision: ---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 1
    res = I; return;
else
    S = varargin{1}; 
end

% determine the interval object
if ~isa(I,'interval')
    temp = I;
    I = S;
    S = temp;
end

% different cases depending on the class of the summand
if isa(S,'interval') && isempty(S)
    % actually holds for all sets, but other checks might be costly
    
    res = interval();
    
elseif isa(S,'interval') || isnumeric(S)

    res = I | S;

elseif isa(S,'zonotope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'polyZonotope') || ...
       isa(S,'mptPolytope') || isa(S,'conPolyZono')

    res = convHull(S,I);

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',I,S));
        
end

%------------- END OF CODE --------------