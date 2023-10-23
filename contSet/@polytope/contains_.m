function res = contains_(P,S,type,tol,varargin)
% contains_ - determines if a polytope contains another set of points
%
% Syntax:
%    res = contains_(P,S,type,tol)
%
% Inputs:
%    P - polytope object
%    S - contSet object, numerical vector
%    tol - numerical tolerance for point in set containment
%
% Outputs:
%    res - true/false whether containment holds true
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]);
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
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper, Viktor Kotsev
% Written:       19-November-2019
% Last update:   26-July-2021 (VG, extended to multiple points)
%                26-April-2022 (added cases for empty objects)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;

% check if P1 is empty
if isemptyobject(P)
    res = false;
    return
end

% check if P2 is empty
if isemptyobject(S)
    res = true;
    return
end

% 1D -> cheap computation of vertices (skip linear program below)
if dim(P) == 1 && ~isnumeric(S)
    res = all(contains_(P,vertices(S),'exact',tol));
    return
end


% get object properties
A = P.A;
b = P.b;

% point in polytope containment
if isnumeric(S)
    
    tmp = A*S - b;
    res = all(tmp < tol | withinTol(tmp,tol));

% other set in polytope containment
else
    % prone to errors, must be carefully checked
    % if isa(S,'polytope') && ~isempty(S.V.val) 
    %     tmp = A*S.V.val - b;
    %     res = all(all(tmp < tol | withinTol(tmp,tol)));
    %     return
    % end

    otherOptions = {};
    if isa(S,'conPolyZono') || isa(S,'polyZonotope')
        otherOptions = {'interval',8,1e-3};
    end

    % loop over all halfspaces
    for i = 1:size(A,1)
        b_ = supportFunc_(S,A(i,:)','upper',otherOptions{:});
        if b_ > b(i) && ~withinTol(b(i),b_,tol)
            return 
        end
    end
    
    res = true;

end

% ------------------------------ END OF CODE ------------------------------
