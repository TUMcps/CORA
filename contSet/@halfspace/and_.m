function res = and_(hs,S,varargin)
% and_ - computes the intersection of a halfspace with a set
%
% Syntax:
%    res = and_(hs,S)
%
% Inputs:
%    hs - halfspace object
%    S - contSet object
%
% Outputs:
%    obj - contSet object
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%    hs = halfspace([1 1],2);
%
%    res = hs & zB;
%
%    figure; hold on;
%    xlim([-1,4]);
%    ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res,[1,2],'FaceColor','g');
%    plot(zB,[1,2],'b','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, conZonotope/and_

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

res = S & hs;

% ------------------------------ END OF CODE ------------------------------
