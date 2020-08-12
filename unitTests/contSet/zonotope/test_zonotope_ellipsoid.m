function res = test_zonotope_ellipsoid
% test_zonotope_rank - unit test function of ellipsoid (and implicitly of
% MVEE, insc_ellipsoid, enc_ellipsoid, MVEE, MVIE)
%
% Syntax:  
%    res = test_zonotope_ellipsoid
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
% See also: insc_ellipsoid,enc_ellipsoid,MVEE,MVIE

% Author:       Victor Gaßmann
% Written:      11-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
%do small dimensions if gurobi is not installed => norm takes much longer
%w/o gurobi
dims = 2:3:4;
dGen = 5;
steps = 2;
rndDirs = 100;
for i=dims
    n = i;
    for j=1:steps
        m = n+j*dGen;
        %randomly generate zonotope
        Z = zonotope([randn(n,1),randn(n,m)]);
        %overapproximations
        Eo_exact = ellipsoid(Z,'o:exact');
        Eo_n = ellipsoid(Z,'o:norm');
        Eo_nb = ellipsoid(Z,'o:norm:bnd');
        %underapproximations
        Eu_exact = ellipsoid(Z,'u:exact');
        Eu_n = ellipsoid(Z,'u:norm');
        %Eu_nb = ellipsoid(Z,'u:norm:bnd');%not implemented yet
        %compute rndDirs random unit directions
        cc = randn(n,rndDirs);
        nn = cc./repmat(sqrt(sum(cc.^2,1)),n,1);
        %check if suppfnc(E)>=suppfnc(Z) (Z in E)
        for k=1:rndDirs
            %overapproximations
            if supportFunc(Eo_exact,nn(:,k))<supportFunc(Z,nn(:,k)) || ...
                    supportFunc(Eo_n,nn(:,k))<supportFunc(Z,nn(:,k)) || ...
               supportFunc(Eo_nb,nn(:,k))<supportFunc(Z,nn(:,k))
                res = false;
                break;
            end
            %underapproximations
            if supportFunc(Eu_exact,nn(:,k))>supportFunc(Z,nn(:,k)) || ...
                    supportFunc(Eu_n,nn(:,k))>supportFunc(Z,nn(:,k))
                res = false;
                break;
            end
        end
        if ~res
            break;
        end
    end
    if ~res
        break;
    end
end
if res
    disp('test_zonotope_ellipsoid successful');
else
    disp('test_zonotope_ellipsoid failed');
end

%------------- END OF CODE --------------
