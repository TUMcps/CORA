function res = isIntersecting(Z,S,varargin)
% isIntersecting - determines if zonotope intersects a set
%
% Syntax:  
%    res = isIntersecting(Z,S)
%    res = isIntersecting(Z,S,type)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope([0 1 1 0;0 1 0 1]);
%    Z2 = zonotope([2 -1 1 0;2 1 0 1]);
%    Z3 = zonotope([3.5 -1 1 0;3 1 0 1]);
% 
%    isIntersecting(Z1,Z2)
%    isIntersecting(Z1,Z3)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_isIntersecting('zonotope',Z,S,varargin{:});

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign values
    Z = vars{1}; S = vars{2};
end


% call function for other set representations
if isa(S,'halfspace') || isa(S,'conHyperplane') || ...
   isa(S,'mptPolytope') || isa(S,'ellipsoid')

    res = isIntersecting(S,Z,varargin{:});

else
    
    res = isIntersecting(conZonotope(Z),S,varargin{:});
    
end

%------------- END OF CODE --------------