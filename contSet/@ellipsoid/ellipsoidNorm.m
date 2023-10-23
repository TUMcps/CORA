function res = ellipsoidNorm(E,p)
% ellipsoidNorm - computes the norm of the point p w.r.t. the
%    ellipsoid-norm induced by the ellipsoid E; this is defined similarly
%    to the zonotope-norm defined in [1, Definition 4].
%
% Syntax:
%    res = ellipsoidNorm(E,p)
%
% Inputs:
%    E - ellipsoid
%    p - nx1-array, with n the dimension of E
%
% Outputs:
%    res - ellipsoid-norm of the point p
%
% Example:
%    q = [0;0];
%    Q = diag([3 2]);
%    E = ellipsoid(Q, q);
%    
%    p = rand([2 1]);
%    
%    d = ellipsoidNorm(E, p);
%
%    figure; hold on
%    plot(E,[1,2],'b');
%    plot(d*E,[1,2],'g');
%    plot(p(1),p(2),'rx');
%
% References:
%    [1] A. Kulmburg, M. Althoff. "On the co-NP-Completeness of the
%        Zonotope Containment Problem", European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       06-July-2021
% Last update:   26-July-2021 (VG, check for degenerate ellipsoid)
%                04-July-2022 (VG, more precise input check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {p,'att','double',{'size',[dim(E),1]}}});

if rank(E) ~= dim(E)
    throw(CORAerror('CORA:degenerateSet',...
        'Not supported for degenerate ellipsoids!'));
end

% The ellipsoid E is given as {x | (x - q)' * Q^(-1) (x - q) <= 1}...
Q = E.Q;
% ... as a consequence, the norm of a point p w.r.t. E would be given by
% p'*Q^(-1)*p, since we have to ignore the center to get a norm.
res = sqrt(p'*(Q\p));
% The square root is just there to make sure that the resulting function is
% a norm (i.e., scaling p by a factor a should yield a|p|, not a^2|p|.

% ------------------------------ END OF CODE ------------------------------
