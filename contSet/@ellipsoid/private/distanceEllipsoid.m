function val = distanceEllipsoid(E1,E2)
% distancePoint - computes the distance from an ellipsoid to an (array of)
% points
%
% Syntax:  
%    D = distancePoint(E,Y)
%
% Inputs:
%    E - ellipsoid object
%    Y - point(s)
%
% Outputs:
%    D - distance(s) between E and Y
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% not supported for degenerate ellipsoids or of different size
if dim(E1) ~= dim(E2)
    error('Dimensions of ellipsoids do not match!');
end
if E1.isdegenerate && E2.isdegenerate
    error('At least one ellipsoid has to be non-degenerate!');
end

max_val = max([svd(E1.Q);svd(E2.Q)]);
E1 = 1/max_val*E1;
E2 = 1/max_val*E2;
q1 = E1.q;
E1 = E1 + (-q1);
E2 = E2 + (-q1);

x2_rem = [];
% make E1 always non-degenerate
if E1.isdegenerate
    tmp = E1;
    E1 = E2;
    E2 = tmp;
end
if E2.isdegenerate
    nt = E2.rank;
    % Q all zero
    if nt==0
        % distance from E1 to center q
        val = max_val*distancePoint(E1,E2.q);
        return;
    end
    [T,~,~] = svd(E2.Q);
    E2 = T'*E2;
    E1 = T'*E1;
    x2_rem = E2.q(nt+1:end);
    % project
    E2 = project(E2,1:nt);
end

n = dim(E1);
nt = dim(E2);

Q1 = inv(E1.Q);
Q2 = inv(E2.Q);
q1 = E1.q;
q2 = E2.q;

x1 = sdpvar(n,1);
x2 = sdpvar(nt,1);

C1 = (x1-q1)'*Q1*(x1-q1) <= 1;
C2 = (x2-q2)'*Q2*(x2-q2) <= 1;

f_obj = norm(x1-[x2;x2_rem]);
diagnostics = optimize([C1;C2],f_obj,sdpsettings('verbose',0));
% no need to backtransform, since distance is invariant under unitary
% transformations
if diagnostics.problem ~= 0
    throw(errOptNumIssue());
end
val = max_val*value(f_obj);
%------------- END OF CODE --------------