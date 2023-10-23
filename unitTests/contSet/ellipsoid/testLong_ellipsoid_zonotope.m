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
% See also: -

% Authors:       Victor Gassmann, Matthias Althoff
% Written:       14-October-2019
% Last update:   07-June-2021 (MA, degenerate case added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
dims = 2:5;
dGen = 5;
steps = 3;
rndDirs = 100;
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
        Zo_n = zonotope(E,m,'outer:norm');
        %Zo_nb = zonotope(E,m,'outer:norm_bnd');
        % degenerate case
        Zo_box_deg = zonotope(E_deg);
        Zo_n_deg = zonotope(E_deg,m,'outer:norm');
        %Zo_nb_deg = zonotope(E_deg,m,'outer:norm_bnd');
        
        %% underapproximations
        % normal case
        Zu_box = zonotope(E,[],'inner:box');
        Zu_n = zonotope(E,m,'inner:norm');
        Zu_nb = zonotope(E,m,'inner:norm_bnd');
        % degenerate case
        Zu_box_deg = zonotope(E_deg,[],'inner:box');
        Zu_n_deg = zonotope(E_deg,m,'inner:norm');
        Zu_nb_deg = zonotope(E_deg,m,'inner:norm_bnd');
        
        
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
            if supportFunc(E,d)>supportFunc(Zo_n,d) || ...
               supportFunc(E,d)>supportFunc(Zo_box,d)
               %supportFunc(E,d)>supportFunc(Zo_nb,d) || ...
                res = false;
                return;
            end
            % degenerate case
            if supportFunc(E_deg,d_deg)>supportFunc(Zo_n_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)>supportFunc(Zo_box_deg,d_deg)
               %supportFunc(E_deg,d_deg)>supportFunc(Zo_nb_deg,d_deg) || ...
                res = false;
                return;
            end
            
            %% underapproximations
            % normal case
            if supportFunc(E,d)<supportFunc(Zu_n,d) || ...
               supportFunc(E,d)<supportFunc(Zu_nb,d) || ...
               supportFunc(E,d)<supportFunc(Zu_box,d)
                res = false;
                return;
            end
            % degenerate case
            if supportFunc(E_deg,d_deg)<supportFunc(Zu_n_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)<supportFunc(Zu_nb_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)<supportFunc(Zu_box_deg,d_deg)
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
