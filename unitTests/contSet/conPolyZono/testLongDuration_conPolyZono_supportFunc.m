function res = testLongDuration_conPolyZono_supportFunc
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
    
    % Random Tests --------------------------------------------------------
    
    % define methods that are tested
    methods = {'split','conZonotope','quadProg','interval'};
    
    % loop over all test cases
    for i = 1:5
        
        % generate random constrained polynomial zonotope
        cPZ = conPolyZono.generateRandom();
        
        % generate random direction
        d = -5 + 10*rand(dim(cPZ),1);
        
        % compute inner-approximation of the bound using linear programming
        cPZ_ = d'*cPZ;
        objFun = @(x) cPZ_.c + sum(cPZ_.G.*prod(x.^cPZ_.expMat,1),2);
        conFun = [];
        if ~isempty(cPZ.A)
            conFun = @(x) deal([],-cPZ_.b + sum(cPZ_.A.* ...
                                  prod(x.^cPZ_.expMat_,1),2));
        end
        p = length(cPZ.id);
        lb = -ones(p,1); ub = ones(p,1);
        options = optimoptions('fmincon','Display','off');
        w = warning(); warning('off');
    
        [~,val_] = fmincon(objFun,zeros(p,1),[],[],[],[],lb,ub,conFun,options);
        
        warning(w);
        
        % loop over all methods
        for j = 1:length(methods)
            
            % compute lower bound
            val = supportFunc(cPZ,d,'lower',methods{j});
            
            % check for correctness
            if val >= val_ + tol
                
                % save variables so that failure can be reproduced
                path = pathFailedTests(mfilename());
                save(path,'cPZ','d');
                
                error('Random test failed!');
            end
        end
    end
end

%------------- END OF CODE --------------