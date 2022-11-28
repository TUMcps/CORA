function res = isIntersecting(hs,S,varargin)
% isIntersecting - determines if a halfspace intersects a set
%
% Syntax:  
%    res = isIntersecting(hs,S)
%    res = isIntersecting(hs,S,type)
%
% Inputs:
%    hs - halfspace object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs = halfspace([-1;-1],0);
%    Z1 = zonotope([0 1 1 0; 0 1 0 1]);
%    Z2 = Z1 - [3;3];
% 
%    isIntersecting(hs,Z1)
%    isIntersecting(hs,Z2)
% 
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(Z1,[1,2],'FaceColor','g');
% 
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(Z2,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/isIntersecting

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  14-Sep-2019
%               20-Nov-2019
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_isIntersecting('halfspace',hs,S,varargin{:});

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign values
    hs = vars{1}; S = vars{2}; type = vars{3};
end


% check user input for correctness
if strcmp(type,'exact')
    if isa(S,'taylm') || isa(S,'polyZonotope') || ... 
        isa(S,'ellipsoid') || isa(S,'capsule') 
        throw(CORAerror('CORA:noops',hs,S));
    end
end

% check for intersection
bound = supportFunc(S,hs.c,'lower');
res = bound <= hs.d;

%------------- END OF CODE --------------