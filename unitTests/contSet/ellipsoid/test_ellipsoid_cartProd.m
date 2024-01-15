function res = test_ellipsoid_cartProd
% test_ellipsoid_cartProd - unit test function of cartProd
%
% Syntax:
%    res = test_ellipsoid_cartProd
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

% Authors:       Victor Gassmann
% Written:       26-July-2021
% Last update:   03-January-2023 (MW, add ellipsoid-numeric case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
load cases.mat E_c 

% empty set
% try 
%     cartProd(E_c{1}.E1,ellipsoid.empty(2));
%     res = false;
% catch ME 
%     if ~strcmp(ME.identifier,'CORA:notSupported')
%         rethrow(ME);
%     else
%         res = true;
%     end
% end

% loop over cases
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    E2 = E_c{i}.E2; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    
    % test non-deg
    Eres_nd = cartProd(E1,E2);
    Y_nd = combineVec(randPoint(E1,2*i),randPoint(E2,2*i));
    if ~contains(Eres_nd,Y_nd)
        res = false;
        break;
    end
    
    % test deg
    Eres_d = cartProd(E1,Ed1);
    Y_d = combineVec(randPoint(E1,2*i),randPoint(Ed1,2*i));
    if ~contains(Eres_d,Y_d)
        res = false;
        break;
    end
    
    % test zero rank ellipsoid
    Eres_0 = cartProd(E1,E0);
    Y_0 = combineVec(randPoint(E1,2*i),randPoint(E0,2*i));
    if ~contains(Eres_0,Y_0)
        res = false;
        break;
    end
end

% ellipsoid-numeric case
E = ellipsoid(eye(2),[1;2]);
num = -1;

% compute Cartesian product
E_ = cartProd(E,num);
E__ = cartProd(num,E);

% check result
if ~compareMatrices(E_.Q,[1 0 0; 0 1 0; 0 0 0]) ...
        || ~compareMatrices(E_.q,[1;2;-1]) ...
        || ~compareMatrices(E__.Q,[0 0 0; 0 1 0; 0 0 1]) ...
        || ~compareMatrices(E__.q,[-1;1;2])
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
