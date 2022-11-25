function memberlist = getMembers(fieldname)
% getMembers - returns the group of admissible values a certain field
%    (param or option) can take
%
% Syntax:
%    memberlist = getMembers(element)
%
% Inputs:
%    fieldname - param or option (char)
%
% Outputs:
%    memberlist - cell array of admissible values for given element
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-January-2021
% Last update:  28-January-2021
% Last revision:---

%------------- BEGIN CODE --------------

if strcmp(fieldname,'R0')
    memberlist = {'capsule','conPolyZono','conZonotope','ellipsoid','interval',...
        'mptPolytope','polyZonotope','probZonotope','zonoBundle','zonotope'};
    
elseif strcmp(fieldname,'U')
    memberlist = {'zonotope','interval'};

elseif strcmp(fieldname,'V')
    memberlist = {'zonotope','ellipsoid'};
    
elseif strcmp(fieldname,'W')
    memberlist = {'zonotope','ellipsoid'};
    
elseif strcmp(fieldname,'Usim')
    memberlist = {'capsule','conPolyZono','conZonotope','ellipsoid','interval',...
        'mptPolytope','polyZonotope','probZonotope','zonoBundle','zonotope'};
    
elseif strcmp(fieldname,'alg')
    memberlist = {'lin','poly','linRem','lin-adaptive','poly-adaptive'};

elseif strcmp(fieldname,'alg4param')
    memberlist = {'lin','poly','linRem'};
    
elseif strcmp(fieldname,'alg4observe')
    memberlist = {'VolMin-A','VolMin-B',...
        'FRad-A','FRad-B','FRad-C',...
        'PRad-A','PRad-B','PRad-C','PRad-D','PRad-E',...
        'CZN-A','CZN-B',...
        'Nom-G','Hinf-G',...
        'ESO-A','ESO-B','ESO-C','ESO-D',...
        'backward','ROPO'};
    
elseif strcmp(fieldname,'algInner')
    memberlist = {'proj','parallelo','scale'};
    
elseif strcmp(fieldname,'reductionTechnique')
    memberlist = {'girard','combastel','pca','methA','methB','methC',...
        'methD','methE','methF','redistribute','cluster','scott','constOpt'};

elseif strcmp(fieldname,'reductionTechnique4nlsys')
    memberlist = [getMembers('reductionTechnique'),...
                  {'approxdep_girard','approxdep_pca'}];
    
elseif strcmp(fieldname,'reductionTechniqueUnderApprox')
    memberlist = {'sum','scale','linProg'};
    
elseif strcmp(fieldname,'linAlg')
    memberlist = {'standard','wrapping-free','fromStart',...
        'decomp','krylov','adaptive'};

elseif strcmp(fieldname,'linAlg4HA')
    memberlist = {'standard','wrapping-free','fromStart','adaptive'};
    
elseif strcmp(fieldname,'guardIntersect')
    memberlist = {'polytope','conZonotope','levelSet',...
        'zonoGirard','pancake','hyperplaneMap','nondetGuard'};
    
elseif strcmp(fieldname,'guardIntersect4enclose')
    % if guardIntersect is one of these, then enclose is mandatory
    memberlist = {'polytope','conZonotope','zonoGirard','nondetGuard'};
    
elseif strcmp(fieldname,'guardIntersect4guardOrder')
    % if guardIntersect is one of these, then guardOrder is mandatory
    memberlist = {'conZonotope','hyperplaneMap'};
    
elseif strcmp(fieldname,'enclose')
    memberlist = {'box','pca','flow'};
    
elseif strcmp(fieldname,'restructureTechnique')
    prefix = {'zonotope','reduceFull','reduce'};
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
    
elseif strcmp(fieldname,'lagrangeRem.simplify')
    memberlist = {'none','simplify','collect','optimize'};
    
elseif strcmp(fieldname,'lagrangeRem.method')
    memberlist = {'interval','taylorModel','zoo'};
 
elseif strcmp(fieldname,'lagrangeRem.zooMethods')
    memberlist = {'interval','affine(int)','affine(bnb)','affine(bnbAdv)',...
        'affine(linQuad)','taylm(int)','taylm(bnb)','taylm(bnbAdv)','taylm(linQuad)'};
    
elseif strcmp(fieldname,'lagrangeRem.optMethod')
    memberlist = {'int','bnb','bnbAdv','linQuad'};
    
elseif strcmp(fieldname,'contractor')
    memberlist = {'linearize','forwardBackward','polyBox'};
    
else
    error("No list of admissible values provided.");
    
end

end

%------------- END OF CODE --------------
