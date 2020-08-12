function res = test_ellipsoid_zonotope
% test_ellipsoid_zonotope - unit test function of zonotope (which
% implicitly is a unit test for insc_zonotope and enc_zonotope)
%
% Syntax:  
%    res = test_ellipsoid_zonotope
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

% Author:       Victor Gaﬂmann
% Written:      14-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
dims = 2:5;
dGen = 5;
steps = 3;
rndDirs = 100;
for i=dims
    %randomly generate ellipsoid
    E = ellipsoid.generateRandom(false,i);
    n = i;
    for j=1:steps
        m = n+j*dGen;
        Zo_box = zonotope(E);
        Zo_n = zonotope(E,m,'o:norm');
        %Zo_nb = zonotope(E,m,'o:norm:bnd');%not implemented yet
        Zu_box = zonotope(E,[],'u:box');
        Zu_n = zonotope(E,m,'u:norm');
        Zu_nb = zonotope(E,m,'u:norm:bnd');
        %compute rndDirs random unit directions
        cc = randn(n,rndDirs);
        nn = cc./repmat(sqrt(sum(cc.^2,1)),n,1);
        %check if suppfnc(E)<=suppfnc(Z) (E in Z)
        for k=1:rndDirs
            d = nn(:,k);
            %overapproximations
            if supportFunc(E,d)>supportFunc(Zo_n,d) || ...
                    supportFunc(E,d)>supportFunc(Zo_box,d)
                res = false;
                break;
            end
            %underapproximations
            if supportFunc(E,d)<supportFunc(Zu_n,d) || ...
                    supportFunc(E,d)<supportFunc(Zu_nb,d) || ...
                    supportFunc(E,d)<supportFunc(Zu_box,d)
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
    disp('test_ellipsoid_zonotope successful');
else
    disp('test_ellipsoid_zonotope failed');
end
%------------- END OF CODE --------------