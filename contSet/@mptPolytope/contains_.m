function res = contains_(P,S,type,tol,varargin)
% contains - determines if a polytope contains another set or a point
%
% Syntax:  
%    res = contains_(P,S)
%    res = contains_(P,S,type)
%    res = contains_(P,S,type,tol)
%
% Inputs:
%    P - mptPolytope object
%    S - contSet object or single point
%    type - always 'exact'
%    tol - numerical tolerance for point in set containment
%
% Outputs:
%    res - true/false
%
% Example: 
%    P1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]);
%    P3 = P2 + [2;0];
%
%    contains(P1,P2)
%    contains(P1,P3)
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P2,[1,2],'g');
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains_

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  26-July-2021 (VG: extended to multiple points)
%               15-November-2022 (MW, return logical array for points)
%               25-November-2022 (MW, rename 'contains')
% Last revision:27-March-2023 (MW, rename contains_)

%------------- BEGIN CODE --------------

% init result
res = false;

% get object properties
A = P.P.A;
b = P.P.b;

% point in polytope containment
if isnumeric(S)
    tmp = A*S - b;
    res = all(tmp < tol | withinTol(tmp,tol));

% other set in polytope containment
else

    % loop over all halfspaces
    for i = 1:size(A,1)
        b_ = supportFunc_(S,A(i,:)','upper','interval',8,1e-3);
        if b_ > b(i) && ~withinTol(b(i),b_,tol)
            return 
        end
    end
    
    res = true;

end

%------------- END OF CODE --------------
