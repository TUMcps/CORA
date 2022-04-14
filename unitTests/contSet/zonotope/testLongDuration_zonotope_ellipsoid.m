function res = testLongDuration_zonotope_ellipsoid
% testLongDuration_zonotope_rank - unit test function of ellipsoid (and implicitly of
% MVEE, insc_ellipsoid, enc_ellipsoid, MVEE, MVIE)
%
% Syntax:  
%    res = testLongDuration_zonotope_ellipsoid
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

% Author:       Victor Gassmann, Matthias Althoff
% Written:      11-October-2019
% Last update:  06-June-2021 (MA, degenerate case added)
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
        %% random zonotopes
        % normal case
        Z = zonotope([randn(n,1),randn(n,m)]);
        % degenerate case
        T = rand(n+1,n);
        Z_deg = T*Z;
        
        %% overapproximations
        % normal case
        Eo_exact = ellipsoid(Z,'o:exact');
        Eo_n = ellipsoid(Z,'o:norm');
        Eo_nb = ellipsoid(Z,'o:norm:bnd');
        % degenerate case
        Eo_exact_deg = ellipsoid(Z_deg,'o:exact');
        Eo_n_deg = ellipsoid(Z_deg,'o:norm');
        Eo_nb_deg = ellipsoid(Z_deg,'o:norm:bnd');
        
        %% underapproximations
        % normal case
        Eu_exact = ellipsoid(Z,'i:exact');
        Eu_n = ellipsoid(Z,'i:norm');
        %Eu_nb = ellipsoid(Z,'i:norm:bnd');
        % degenerate case
        Eu_exact_deg = ellipsoid(Z_deg,'i:exact');
        Eu_n_deg = ellipsoid(Z_deg,'i:norm');
        %Eu_nb_deg = ellipsoid(Z_deg,'i:norm:bnd');
        
        %% compute rndDirs random unit directions
        % normal case
        cc = randn(n,rndDirs);
        nn = cc./repmat(sqrt(sum(cc.^2,1)),n,1);
        % degenerate case
        nn_deg = T*nn;
        
        % check if suppfnc(E)>=suppfnc(Z) (Z in E)
        for k=1:rndDirs
            %% overapproximations
            % normal case
            if supportFunc(Eo_exact,nn(:,k))<supportFunc(Z,nn(:,k)) || ...
               supportFunc(Eo_n,nn(:,k))<supportFunc(Z,nn(:,k)) || ...
               supportFunc(Eo_nb,nn(:,k))<supportFunc(Z,nn(:,k))
                res = false;
                break;
            end
            % degenerate case
            if supportFunc(Eo_exact_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eo_n_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eo_nb_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k))
                res = false;
                break;
            end
            
            %% underapproximations
            % normal case
            if supportFunc(Eu_exact,nn(:,k))>supportFunc(Z,nn(:,k)) || ...
               supportFunc(Eu_n,nn(:,k))>supportFunc(Z,nn(:,k)) %|| ...
               %supportFunc(Eu_nb,nn(:,k))>supportFunc(Z,nn(:,k)) 
                res = false;
                break;
            end
            % degenerate case
            if supportFunc(Eu_exact_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eu_n_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) %|| ...
               %supportFunc(Eu_nb_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) 
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
    disp('testLongDuration_zonotope_ellipsoid successful');
else
    disp('testLongDuration_zonotope_ellipsoid failed');
end

%------------- END OF CODE --------------
