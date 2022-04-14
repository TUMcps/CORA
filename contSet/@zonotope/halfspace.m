function Z = halfspace(Z)
% halfspace - generates the halfspace representation of the zonotope,
%    which is stored in the zonotope object (Z.halfspace)
%
% Syntax:  
%    Z = halfspace(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object, including halfspace representation
%
% Example: 
%    Z = zonotope([0;0],rand(2,4));
%    Z = halfspace(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  06-April-2017
%               16-October-2019: added box case (VG)
% Last revision:---

%------------- BEGIN CODE --------------

%check if zonotope is a parallelepiped
G = generators(Z);
if all(size(G)==size(G,1)) && rank(G)==size(G,1)
    c = center(Z);
    n = size(G,1);
    %x=c+G*u, |u|<=1 -> [I;-I]*inv(G)*(x-c)<=1
    A = [eye(n);-eye(n)];
    H = A*inv(G);
    K = ones(2*n,1)+H*c;
    %normalize H to H(i,:)=1
    hn = sqrt(sum(H.^2,2));
    H = H./repmat(hn,1,size(H,2));
    K = K./hn;
else
    %convert zonotope to polytope and retrieve halfspace representation
    P = polytope(Z);
    H = get(P,'H');
    K = get(P,'K');
end
%write to object structure
Z.halfspace.H=H;
Z.halfspace.K=K;
Z.halfspace.equations=length(K);

%------------- END OF CODE --------------