function res = test_conPolyZono_supportFunc
% test_conPolyZono_supportFunc - unit test function for support function 
%                                enclosure of constrained polynomial 
%                                zonotopes
%
% Syntax:  
%    test_conPolyZono_supportFunc
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/supportFunc

% Author:       Niklas Kochdumper
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
	
    res = 1;
    tol = 1e-6;

    % Analytical Tests ----------------------------------------------------

    % define constrained polynomial zonotope
    c = [0;0];
    G = [1 0 1 -1; 0 2 1 2];
    expMat = [1 2 1 0; 0 0 1 2; 0 0 0 0];
    A = [1 1 0.5];
    b = 0.5;
    expMat_ = [0 1 0;1 0 0; 0 0 1];

    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
    
    % define direction
    d = [1;1];

    % remove constraints
    val = supportFunc(cPZ,d,'lower','quadProg');
    
    % compute exact solution
    temp = d'*G;
    H = 2*blkdiag([temp(2) 0.5*temp(3);0.5*temp(3) temp(4)],0);
    f = [temp(1);0;0];
    Aeq = [A(2);A(1);A(3)];
    
    options = optimoptions('quadprog','Display','off');
    [~,val_] = quadprog(H,f,[],[],Aeq',b,-ones(3,1),ones(3,1),[],options);
    
%     % visualization
%     figure; hold on
%     plot(cPZ,[1,2],'r','Splits',15);
%     ch = conHyperplane(d,val_);
%     plot(ch,[1,2],'b');

    % compare with exact solution
    if abs(val-val_) >= tol
        error('Analytical test failed!');
    end
    
end

%------------- END OF CODE --------------