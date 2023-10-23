function res = isIntersecting_(hs,S,type,varargin)
% isIntersecting_ - determines if a halfspace intersects a set
%
% Syntax:
%    res = isIntersecting_(hs,S)
%    res = isIntersecting_(hs,S,type)
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
% See also: contSet/isIntersecting, conHyperplane/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       16-May-2018
% Last update:   14-September-2019
%                20-November-2019
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% halsfspace forms a special case
if isa(S,'halfspace')
    res = sum(abs(S.c/norm(S.c) + hs.c/norm(hs.c))) > 1e-10 | ...
                                S.d/norm(S.c) + hs.d/norm(hs.c) > 0;
    return;
end

% check user input for correctness
if strcmp(type,'exact')
    if isa(S,'taylm') || isa(S,'polyZonotope') || ... 
        isa(S,'ellipsoid') || isa(S,'capsule') 
        throw(CORAerror('CORA:noops',hs,S));
    end
end

% check for intersection
bound = supportFunc_(S,hs.c,'lower','interval',8,1e-3);
res = bound <= hs.d;

% ------------------------------ END OF CODE ------------------------------
