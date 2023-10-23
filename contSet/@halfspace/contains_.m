function res = contains_(hs,S,varargin)
% contains_ - determines if a halfspace contains a set or a point
%
% Syntax:
%    res = contains_(hs,S)
%
% Inputs:
%    hs - halfspace object
%    S - contSet object or single point
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs = halfspace([-1;-1],0);
%    Z1 = zonotope([3 1 1 0;3 1 0 1]);
%    Z2 = Z1 - [3;3];
%
%    contains(hs,Z1)
%    contains(hs,Z2)
%
%    figure; hold on;
%    xlim([-6,6]); ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(Z1,[1,2],'FaceColor','g');
%
%    figure; hold on;
%    xlim([-6,6]); ylim([-6,6]);
%    plot(hs,[1,2],'b');
%    plot(Z2,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       14-June-2016
% Last update:   27-July-2016
%                02-September-2019
%                19-November-2019 (NK, extend to all set representations)
%                25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% point in halfspace containment
if isnumeric(S)
    tmp = hs.c' * S;
    res = tmp < hs.d | withinTol(tmp,hs.d);
    
% set in halfspace containment
else
    val = supportFunc_(S,hs.c,'upper','interval',8,1e-3);
    res = val < hs.d | withinTol(val,hs.d);
end

% ------------------------------ END OF CODE ------------------------------
