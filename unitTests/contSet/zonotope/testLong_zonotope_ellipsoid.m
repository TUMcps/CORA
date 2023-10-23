function res = testLong_zonotope_ellipsoid
% testLong_zonotope_ellipsoid - unit test function of ellipsoid
%    (and implicitly of MVEE, insc_ellipsoid, enc_ellipsoid, MVEE, MVIE)
%
% Syntax:
%    res = testLong_zonotope_ellipsoid
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: insc_ellipsoid,enc_ellipsoid,MVEE,MVIE

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       11-October-2019
% Last update:   06-June-2021 (MA, degenerate case added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
        Eo_exact = ellipsoid(Z,'outer:exact');
        Eo_n = ellipsoid(Z,'outer:norm');
        Eo_nb = ellipsoid(Z,'outer:norm_bnd');
        % degenerate case
        Eo_exact_deg = ellipsoid(Z_deg,'outer:exact');
        Eo_n_deg = ellipsoid(Z_deg,'outer:norm');
        Eo_nb_deg = ellipsoid(Z_deg,'outer:norm_bnd');
        
        %% underapproximations
        % normal case
        Eu_exact = ellipsoid(Z,'inner:exact');
        Eu_n = ellipsoid(Z,'inner:norm');
        %Eu_nb = ellipsoid(Z,'inner:norm_bnd');
        % degenerate case
        Eu_exact_deg = ellipsoid(Z_deg,'inner:exact');
        Eu_n_deg = ellipsoid(Z_deg,'inner:norm');
        %Eu_nb_deg = ellipsoid(Z_deg,'inner:norm_bnd');
        
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
                throw(CORAerror('CORA:testFailed'));
            end
            % degenerate case
            if supportFunc(Eo_exact_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eo_n_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eo_nb_deg,nn_deg(:,k))<supportFunc(Z_deg,nn_deg(:,k))
                throw(CORAerror('CORA:testFailed'));
            end
            
            %% underapproximations
            % normal case
            if supportFunc(Eu_exact,nn(:,k))>supportFunc(Z,nn(:,k)) || ...
               supportFunc(Eu_n,nn(:,k))>supportFunc(Z,nn(:,k)) %|| ...
               %supportFunc(Eu_nb,nn(:,k))>supportFunc(Z,nn(:,k)) 
                throw(CORAerror('CORA:testFailed'));
            end
            % degenerate case
            if supportFunc(Eu_exact_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) || ...
               supportFunc(Eu_n_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) %|| ...
               %supportFunc(Eu_nb_deg,nn_deg(:,k))>supportFunc(Z_deg,nn_deg(:,k)) 
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
