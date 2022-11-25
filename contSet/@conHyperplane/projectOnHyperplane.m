function res = projectOnHyperplane(h, S)
% projectOnHyperplane - projects a set onto a hyperplane
%
% Syntax:  
%    S = projectOnHyperplane(h, S)
%
% Inputs:
%    h - conHyperplane object
%    S - contSet object
%
% Outputs:
%    res - projected set
%
% Example: 
%    h = conHyperplane([1 1],1);
%    zono = zonotope([2 1 -1;2 0 1]);
%
%    res = projectOnHyperplane(h,zono);
%
%    figure
%    hold on
%    xlim([-3,5]);
%    ylim([-3,5]);
%    plot(h,[1,2],'r');
%    plot(zono,[1,2],'g');
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Author:       Niklas Kochdumper
% Written:      13-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get object properties
    c = h.h.c;
    d = h.h.d;

    % normalization
    temp = norm(c);

    c = c./temp;
    d = d/temp;

    % linear map A*x + b for the projection
    A = eye(length(c)) - c*c';
    b = d*c;

    % project the set
    res = A*S + b;

end