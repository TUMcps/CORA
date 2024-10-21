function res = testLong_ellipsoid_zonotope
% testLong_ellipsoid_zonotope - unit test function of zonotope
%    (which implicitly is a unit test for insc_zonotope and enc_zonotope)
%
% Syntax:
%    res = testLong_ellipsoid_zonotope
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
% See also: none

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       14-October-2019
% Last update:   07-June-2021 (MA, degenerate case added)
%                22-April-2024 (TL, added tolerance)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
dims = 2:5;
dGen = 5;
steps = 3;
rndDirs = 100;
TOL = 1e-6;

for i=dims
    
    %%% generate all variables necessary to replicate results
    n = i;
    %% random ellipsoids
    % normal case
    E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
    % degenerate case
    T = rand(n+1,n);
    %compute rndDirs random directions
    cc = randn(n,rndDirs);
    %%%
    E_deg = T*E;
  
    for j=1:steps
        m = n+j*dGen;
        
        %% overapproximations
        % normal case
        Zo_box = zonotope(E);
        Zo_n = zonotope(E,'outer:norm',m);
        %Zo_nb = zonotope(E,'outer:norm_bnd',m);
        % degenerate case
        Zo_box_deg = zonotope(E_deg);
        Zo_n_deg = zonotope(E_deg,'outer:norm',m);
        %Zo_nb_deg = zonotope(E_deg,'outer:norm_bnd',m);
        
        %% underapproximations
        % normal case
        Zu_box = zonotope(E,'inner:box');
        Zu_n = zonotope(E,'inner:norm',m);
        Zu_nb = zonotope(E,'inner:norm_bnd',m);
        % degenerate case
        Zu_box_deg = zonotope(E_deg,'inner:box');
        Zu_n_deg = zonotope(E_deg,'inner:norm',m);
        Zu_nb_deg = zonotope(E_deg,'inner:norm_bnd',m);
        
        
        nn = cc./repmat(sqrt(sum(cc.^2,1)),n,1);
        %check if suppfnc(E)<=suppfnc(Z) (E in Z)
        for k=1:rndDirs
            
            %% directions
            % normal case
            d = nn(:,k);
            % degenerate case
            d_deg = T*nn(:,k);
            
            %% overapproximations
            % normal case
            assertLoop(supportFunc(E,d) <= supportFunc(Zo_n,d)+TOL,i,j,k)
            assertLoop(supportFunc(E,d) <= supportFunc(Zo_box,d)+TOL,i,j,k)
            % assertLoop(supportFunc(E,d) <= supportFunc(Zo_nb,d)+TOL,i,j,k)

            % degenerate case
            assertLoop(supportFunc(E_deg,d_deg) <= supportFunc(Zo_n_deg,d_deg)+TOL,i,j,k)
            assertLoop(supportFunc(E_deg,d_deg) <= supportFunc(Zo_box_deg,d_deg)+TOL,i,j,k)
            %assertLoop(supportFunc(E_deg,d_deg) <= supportFunc(Zo_nb_deg,d_deg)+TOL,i,j,k)
            
            %% underapproximations
            % normal case
            assertLoop(supportFunc(E,d)+TOL >= supportFunc(Zu_n,d),i,j,k)
            assertLoop(supportFunc(E,d)+TOL >= supportFunc(Zu_nb,d),i,j,k)
            assertLoop(supportFunc(E,d)+TOL >= supportFunc(Zu_box,d),i,j,k)

            % degenerate case
            assertLoop(supportFunc(E_deg,d_deg)+TOL >= supportFunc(Zu_n_deg,d_deg),i,j,k)
            assertLoop(supportFunc(E_deg,d_deg)+TOL >= supportFunc(Zu_nb_deg,d_deg),i,j,k)
            assertLoop(supportFunc(E_deg,d_deg)+TOL >= supportFunc(Zu_box_deg,d_deg),i,j,k)
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
