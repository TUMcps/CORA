function res = testLongDuration_ellipsoid_zonotope
% testLongDuration_ellipsoid_zonotope - unit test function of zonotope
%    (which implicitly is a unit test for insc_zonotope and enc_zonotope)
%
% Syntax:  
%    res = testLongDuration_ellipsoid_zonotope
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

% Author:       Victor Gassmann, Matthias Althoff
% Written:      14-October-2019
% Last update:  07-June-2021 (MA, degenerate case added)
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
dims = 2:5;
dGen = 5;
steps = 3;
rndDirs = 100;
for i=dims
    n = i;
    
    %% random ellipsoids
    % normal case
    E = ellipsoid.generateRandom(false,i);
    % degenerate case
    T = rand(n+1,n);
    E_deg = T*E;
  
    for j=1:steps
        m = n+j*dGen;
        
        %% overapproximations
        % normal case
        Zo_box = zonotope(E);
        Zo_n = zonotope(E,m,'o:norm');
        %Zo_nb = zonotope(E,m,'o:norm:bnd');
        % degenerate case
        Zo_box_deg = zonotope(E_deg);
        Zo_n_deg = zonotope(E_deg,m,'o:norm');
        %Zo_nb_deg = zonotope(E_deg,m,'o:norm:bnd');
        
        %% underapproximations
        % normal case
        Zu_box = zonotope(E,[],'i:box');
        Zu_n = zonotope(E,m,'i:norm');
        Zu_nb = zonotope(E,m,'i:norm:bnd');
        % degenerate case
        Zu_box_deg = zonotope(E_deg,[],'i:box');
        Zu_n_deg = zonotope(E_deg,m,'i:norm');
        Zu_nb_deg = zonotope(E_deg,m,'i:norm:bnd');
        
        
        %compute rndDirs random unit directions
        cc = randn(n,rndDirs);
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
                break;
            end
            % degenerate case
            if supportFunc(E_deg,d_deg)>supportFunc(Zo_n_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)>supportFunc(Zo_box_deg,d_deg)
               %supportFunc(E_deg,d_deg)>supportFunc(Zo_nb_deg,d_deg) || ...
                res = false;
                break;
            end
            
            %% underapproximations
            % normal case
            if supportFunc(E,d)<supportFunc(Zu_n,d) || ...
               supportFunc(E,d)<supportFunc(Zu_nb,d) || ...
               supportFunc(E,d)<supportFunc(Zu_box,d)
                res = false;
                break;
            end
            % degenerate case
            if supportFunc(E_deg,d_deg)<supportFunc(Zu_n_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)<supportFunc(Zu_nb_deg,d_deg) || ...
               supportFunc(E_deg,d_deg)<supportFunc(Zu_box_deg,d_deg)
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
    disp('testLongDuration_ellipsoid_zonotope successful');
else
    disp('testLongDuration_ellipsoid_zonotope failed');
end
%------------- END OF CODE --------------