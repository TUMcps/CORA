function E = convHull(E1,E2)
% convHull - computes and over-approximation of the convex hull of two
%            ellipsoids
%
% Syntax:  
%    E = enclose(E1,E2)
%
% Inputs:
%    E1 - first ellipsoid object
%    E2 - second ellipsoid object
%
% Outputs:
%    E - ellipsoid that encloses the convex hull of E1 and E2
%
% Example: 
%    E1 = ellipsoid(eye(2));
%    E2 = ellipsoid([1,0;0,3],[1;-1]);
%
%    E = convHull(E1,E2);
%
%    figure
%    hold on
%    plot(E1,[1,2],'r');
%    plot(E2,[1,2],'b');
%    plot(E,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check input arguments
if nargin~=2 || ( ~isa(E1,'ellipsoid') || ~isa(E2,'ellipsoid'))
    error('Wrong type of arguments');
end
if length(E1.Q)~=length(E2.Q)
    error('Ellipsoids have to have same dimensions');
end
if E1.dim<length(E1.Q) && E2.dim<length(E2.Q)
    warning('Two degenerate ellipsoids: We assume that the convex hull is full-dimensional');
end
if length(E1.Q)==1
    error('not implemented for one-dimensional ellipsoids');
end
q1 = E1.q;
q2 = E2.q;
Q1 = E1.Q;
Q2 = E2.Q;
n = length(Q1);
% %Taken from https://www2.isye.gatech.edu/~nemirovs/Lect_ModConvOpt.pdf,
% %3.7.3.2
[V1,D1] = eig(Q1);
[V2,D2] = eig(Q2);
A1 = V1*sqrt(D1);
A2 = V2*sqrt(D2);
Y = sdpvar(n,n);%symmetric assumption is fine here
z = sdpvar(n,1);
lbda1 = sdpvar;
lbda2 = sdpvar;
t = sdpvar;

M1 = [eye(n),Y*q1-z,Y*A1;
     (Y*q1-z)',1-lbda1,zeros(1,n);
     A1'*Y,zeros(n,1),lbda1*eye(n)];
M2 = [eye(n),Y*q2-z,Y*A2;
     (Y*q2-z)',1-lbda2,zeros(1,n);
     A2'*Y,zeros(n,1),lbda2*eye(n)];
Cnts = [t<=geomean(Y),M1>=0,M2>=0,Y>=0];
optimize(Cnts,-t,sdpsettings('verbose',0));

Q = inv(value(Y)^2);
q = value(Y)\value(z);
E = ellipsoid(Q,q);

%------------- END OF CODE --------------