function res = contains(hs,S)
% contains - determines if a halfspace contains a set or a point
%
% Syntax:  
%    res = contains(hs,S)
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
% See also: zonotope/contains

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-June-2016
% Last update:  27-July-2016
%               02-September-2019
%               19-November-2019 (NK, extend to all set representations)
%               25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_contains('halfspace',hs,S);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    hs = vars{1}; S = vars{2};
end


% point in halfspace containment
if isnumeric(S)
    tmp = hs.c' * S;
    res = tmp < hs.d | withinTol(tmp,hs.d);
    
% set in halfspace containment
else
    val = supportFunc(S,hs.c,'upper');
    res = val < hs.d | withinTol(val,hs.d);
end

%------------- END OF CODE --------------