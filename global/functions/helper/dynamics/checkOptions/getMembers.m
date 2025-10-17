function memberlist = getMembers(fieldname)
% getMembers - returns the group of admissible values a certain field
%    (param or option) can take
%
% Syntax:
%    memberlist = getMembers(element)
%
% Inputs:
%    fieldname - field name in params or options struct (char)
%
% Outputs:
%    memberlist - cell array of admissible values for given element
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Matthias Althoff
% Written:       26-January-2021
% Last update:   28-January-2021
%                05-January-2022 (MA, type constrained simulation added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch fieldname
    case 'R0'
        memberlist = {'capsule','conPolyZono','conZonotope','ellipsoid','interval',...
            'polytope','polyZonotope','probZonotope','spectraShadow','zonoBundle','zonotope'};
    
    case 'R0conf'
        memberlist = {'zonotope'};

    case 'Rend'
        memberlist = {'polytope','interval'};
    
        % disturbances
    case 'U'
        memberlist = {'zonotope','interval','ellipsoid','polyZonotope',...
            'conZonotope','capsule','polytope','conPolyZono','zonoBundle'};

    case 'Usim'
        memberlist = {'capsule','conPolyZono','conZonotope','ellipsoid','interval',...
            'polytope','polyZonotope','probZonotope','zonoBundle','zonotope'};

    case 'V'
        memberlist = {'zonotope','ellipsoid','interval'};
    
    case 'W'
        memberlist = {'zonotope','ellipsoid','interval'};

        % specification
    case 'safeSet'
        memberlist = {'polytope'};

    case 'unsafeSet'
        memberlist = {'polytope'};

        % alg
    case 'alg'
        memberlist = {'lin','poly','linRem','lin-adaptive','poly-adaptive'};

    case 'alg4DT' % discrete time
        memberlist = {'lin','lin-adaptive','poly-adaptive'};
    
    case 'alg4DA'
        memberlist = {'lin','lin-adaptive'};
    
    case 'alg4param'
        memberlist = {'lin','poly','linRem'};
    
    case 'alg4observe'
        memberlist = {'VolMin-A','VolMin-B',...
            'FRad-A','FRad-B','FRad-C',...
            'PRad-A','PRad-B','PRad-C','PRad-D','PRad-E',...
            'CZN-A','CZN-B',...
            'Nom-G','Hinf-G',...
            'ESO-A','ESO-B','ESO-C','ESO-D',...
            'backward','ROPO'};
        
    case 'algInner' % inner
        memberlist = {'proj','parallelo','scale','minkdiff'};
        
    case {'linAlg','reachAlg'} % allow same algorithms
        memberlist = {'standard','wrapping-free','fromStart',...
            'decomp','krylov','adaptive','supportFunc'};
   
    case 'linAlg4backward'
        memberlist = {'inner:EA:timepoint','outer:EA:timepoint','inner:EA:timeinterval',...
                      'inner:AE:timepoint','outer:AE:timepoint','outer:AE:timeinterval'};
        %%% memberlist = {'minimal','maximal','minimal:naive'};
    
    case 'linAlg4HA'
        memberlist = {'standard','wrapping-free','fromStart','adaptive'};

        % reduction
    case 'reductionTechnique'
        memberlist = {'girard','combastel','pca','methA','methB','methC',...
            'methD','methE','methF','redistribute','cluster','scott','constOpt'};
    
    case 'reductionTechnique4nlsys'
        memberlist = [getMembers('reductionTechnique'),...
                      {'approxdep_girard','approxdep_pca'}];
        
    case 'reductionTechniqueUnderApprox'
        memberlist = {'sum','scale','linProg'};
        
        % guard
    case 'guardIntersect'
        memberlist = {'polytope','conZonotope','levelSet',...
            'zonoGirard','pancake','hyperplaneMap','nondetGuard'};
        
    case 'guardIntersect4enclose'
        % if guardIntersect is one of these, then enclose is mandatory
        memberlist = {'polytope','conZonotope','zonoGirard','nondetGuard'};
        
    case 'guardIntersect4guardOrder'
        % if guardIntersect is one of these, then guardOrder is mandatory
        memberlist = {'conZonotope','hyperplaneMap'};
        
    case 'enclose'
        memberlist = {'box','pca','flow'};
        
    case 'restructureTechnique'
        prefix = {'zonotope','reduceFull','reducePart','reduceDI','reduce'};
        redMethod = getMembers('reductionTechnique');
        memberlist = {};
        for i=1:length(redMethod)
            redMethod{i}(1) = upper(redMethod{i}(1));
        end
        for i=1:length(prefix)
            for j=1:length(redMethod)
                memberlist = [memberlist, [prefix{i} redMethod{j}]];
            end
        end
        
        % lagrangeRem
    case 'lagrangeRem.simplify'
        memberlist = {'none','simplify','collect','optimize'};
        
    case 'lagrangeRem.method'
        memberlist = {'interval','taylorModel','zoo'};
     
    case 'lagrangeRem.zooMethods'
        memberlist = {'interval','affine(int)','affine(bnb)','affine(bnbAdv)',...
            'affine(linQuad)','taylm(int)','taylm(bnb)','taylm(bnbAdv)','taylm(linQuad)'};
        
    case 'lagrangeRem.optMethod'
        memberlist = {'int','bnb','bnbAdv','linQuad'};

        % ---
    case 'contractor'
        memberlist = {'linearize','forwardBackward','polyBox'};
        
    case 'type'
        memberlist = {'standard','gaussian','rrt','constrained'};
        
    case 'norm'
        memberlist = {'interval','frob'};
    
    case 'armaxAlg'
        memberlist = {'exactAddition', 'tvpEfficient', 'tvpGeneral'};   
    
    case 'Y0'
        memberlist = {'zonotope','interval','ellipsoid','polyZonotope',...
            'conZonotope','capsule','polytope','conPolyZono','zonoBundle'}; 
        
    case 'idAlg'
        memberlist = {'gp', 'cgp'};   

        % cs
    case 'cs.cost'
        memberlist = {'interval','frob'};
    
    case 'cs.constraints'
        memberlist = {'half','gen'};

    case 'cs.task'
        memberlist = {'reg','class'};

    case 'cs.outMethod'
        memberlist = {'','RMSE', 'search', 'searchG', 'MILP'};

    case 'cs.recMethod'
        memberlist = {'','recUnSiCo', 'recAlpha'};

    otherwise
        throw(CORAerror('CORA:wrongValue','first','Check file.'));

end

end

% ------------------------------ END OF CODE ------------------------------
