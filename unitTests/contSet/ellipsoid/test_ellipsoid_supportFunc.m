function res = test_ellipsoid_supportFunc
% test_ellipsoid_supportFunc - unit test function of supportFunc
%
% Syntax:  
%    res = test_ellipsoid_supportFunc
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gaßmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res=true;
E = ellipsoid.generateRandom();
B = boundary(E,1000);
[V,D] = eig(E.Q);
[~,ind] = sort(diag(D),'descend');
D = D(ind,ind);
V = V(:,ind);
if ~E.isdegenerate
    for i=1:length(B)
        b = B(:,i);
        if (1-(b-E.q)'*inv(E.Q)*(b-E.q))>E.TOL
            res = false;
            break;
        end
    end
else
    Vp = V(:,1:E.dim);
    Vn = V(:,E.dim+1:end);
    Qp = D(1:E.dim,1:E.dim);
    for i=1:length(B)
        b = B(:,i);
        b0 = b-E.q;
        if ~all(abs(Vn'*b0)<=E.TOL)
            res = false;
            break;
        end
        bp = Vp'*b0;
        if abs(1-bp'*inv(Qp)*bp)>E.TOL
            res = false;
            break;
        end
    end
end
if res
    disp('test_ellipsoid_supportFunc successful');
else
    disp('test_ellipsoid_supportFunc failed');
end

%------------- END OF CODE --------------
