function condfun = getCondfunDynParameterOptions(field,listname)
% getCondfunDynParameterOptions - get the condition function of dynamic parameters
%
% Syntax:
%    res = getCondfunDynParameterOptions(field)
%
% Inputs:
%    field - struct field in options
%
% Outputs:
%    condfun - condition function for the dynamic parameters, or empty
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getCondfunDynParameter

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for condition functions in in options
switch field
    case 'error'
        condfun = @aux_getCondfunOptions_error;
    case 'timeStep'
        condfun = @aux_getCondfunOptions_timeStep;
    case 'taylorTerms'
        condfun = @aux_getCondfunOptions_taylorTerms;
    case 'zonotopeOrder'
        condfun = @aux_getCondfunOptions_zonotopeOrder;
    case 'tensorOrder'
        condfun = @aux_getCondfunOptions_tensorOrder;
    case 'l'
        condfun = @aux_getCondfunOptions_l;
    case 'partition'
        condfun = @aux_getCondfunOptions_partition;
    case 'krylovError'
        condfun = @aux_getCondfunOptions_krylovError;
    case 'krylovStep'
        condfun = @aux_getCondfunOptions_krylovStep;
    case 'errorOrder'
        condfun = @aux_getCondfunOptions_errorOrder;
    case 'errorOrder3'
        condfun = @aux_getCondfunOptions_errorOrder3;
    case 'vertSamp'
        condfun = @aux_getCondfunOptions_vertSamp;
    case 'stretchFac'
        condfun = @aux_getCondfunOptions_stretchFac;
    case 'R'
        condfun = @aux_getCondfunOptions_R;
    case 'nrConstInp'
        condfun = @aux_getCondfunOptions_nrConstInp;
    case 'fracInpVert'
        condfun = @aux_getCondfunOptions_fracInpVert;
    case 'fracVert'
        condfun = @aux_getCondfunOptions_fracVert;
    case 'p_conf'
        condfun = @aux_getCondfunOptions_p_conf;
    case 'enclose'
        condfun = @aux_getCondfunOptions_enclose;
    case 'guardOrder'
        condfun = @aux_getCondfunOptions_guardOrder;
    case 'lagrangeRem.zooMethods'
        condfun = @aux_getCondfunOptions_lagrangeRem_zooMethods;
    case 'lagrangeRem.optMethod'
        condfun = @aux_getCondfunOptions_lagrangeRem_optMethod;
    case 'lagrangeRem.maxOrder'
        condfun = @aux_getCondfunOptions_lagrangeRem_maxOrder;
    case 'lagrangeRem.tolerance'
        condfun = @aux_getCondfunOptions_lagrangeRem_tolerance;
    case 'lagrangeRem.eps'
        condfun = @aux_getCondfunOptions_lagrangeRem_eps;
    case 'intermediateOrder'
        condfun = @aux_getCondfunOptions_intermediateOrder;
    case 'polyZono.maxDepGenOrder'
        condfun = @aux_getCondfunOptions_polyZono_maxDepGenOrder;
    case 'polyZono.maxPolyZonoRatio'
        condfun = @aux_getCondfunOptions_polyZono_maxPolyZonoRatio;
    case 'polyZono.restructureTechnique'
        condfun = @aux_getCondfunOptions_polyZono_restructureTechnique;
    case 'inpChanges'
        condfun = @aux_getCondfunOptions_inpChanges;

    otherwise
        % no condition is given
        condfun = [];
end
    
end


% Auxiliary functions -----------------------------------------------------

% options.<field> ---------------------------------------------------------

% error
function res = aux_getCondfunOptions_error(sys,func,params,options)
    res = strcmp(options.linAlg,'adaptive');
end

% timeStep
function res = aux_getCondfunOptions_timeStep(sys,func,params,options)
    if isa(sys,'linParamSys')
        res = false;
    elseif isfield(options, 'linAlg')
        res = ~strcmp(options.linAlg,'adaptive');
    elseif isfield(options, 'alg')
        res = ~contains(options.alg,'adaptive');
    else
        % timeStep is required
        res = true;
    end
end

% taylorTerms
function res = aux_getCondfunOptions_taylorTerms(sys,func,params,options)
    if isa(sys,'linParamSys')
        res = false;
    elseif isfield(options, 'linAlg')
        res = ~strcmp(options.linAlg,'adaptive');
    elseif isfield(options, 'alg')
        res = ~contains(options.alg,'adaptive');
    else
        % taylorTerms is required
        res = true;
    end
end

% zonotopeOrder
function res = aux_getCondfunOptions_zonotopeOrder(sys,func,params,options)
    if isa(sys,'linParamSys')
        res = false;
    elseif isfield(options, 'linAlg')
        res = ~strcmp(options.linAlg,'adaptive');
    elseif isfield(options, 'alg')
        res = ~contains(options.alg,'adaptive');
    else
        % zonotopeOrder is required
        res = true;
    end
end

% tensorOrder
function res = aux_getCondfunOptions_tensorOrder(sys,func,params,options)
    if startsWith(class(sys),'lin')
        res = false;
    elseif isfield(options,'alg')
        res = ~contains(options.alg,'adaptive');
    else
        res = true;
    end
end

% l
function res = aux_getCondfunOptions_l(sys,func,params,options)
    res = strcmp(options.linAlg,'supportFunc');
end

% partition
function res = aux_getCondfunOptions_partition(sys,func,params,options)
    res = strcmp(options.linAlg,'decomp');
end

% krylovError
function res = aux_getCondfunOptions_krylovError(sys,func,params,options)
    res = strcmp(options.linAlg,'krylov');
end

% krylovStep
function res = aux_getCondfunOptions_krylovStep(sys,func,params,options)
    res = strcmp(options.linAlg,'krylov');
end

% errorOrder
function res = aux_getCondfunOptions_errorOrder(sys,func,params,options)
    if strcmp(func,'reachInnerScaling')
        res = true;
    elseif strcmp(func,'reachInnerParallelotope')
        res = options.tensorOrder > 2;
    else
        res = ~startsWith(class(sys),'lin') && ~contains(options.alg,'adaptive') && options.tensorOrder>2;
    end
end

% errorOrder3
function res = aux_getCondfunOptions_errorOrder3(sys,func,params,options)
    res = ~contains(options.alg,'adaptive') && options.tensorOrder>3;
end

% nrConstInp
function res = aux_getCondfunOptions_nrConstInp(sys,func,params,options)
    res = any(strcmp(options.type,{'standard','gaussian'}));
end

% fracInpVert
function res = aux_getCondfunOptions_fracInpVert(sys,func,params,options)
    res = any(strcmp(options.type,{'standard'}));
end

% fracVert
function res = aux_getCondfunOptions_fracVert(sys,func,params,options)
    res = any(strcmp(options.type,{'standard','constrained'}));
end

% p_conf
function res = aux_getCondfunOptions_p_conf(sys,func,params,options)
    res = ~strcmp(func,'simulateRandom') || any(strcmp(options.type,{'gaussian'}));
end

% vertSamp
function res = aux_getCondfunOptions_vertSamp(sys,func,params,options)
    res = strcmp(options.type,'rrt');
end

% stretchFac
function res = aux_getCondfunOptions_stretchFac(sys,func,params,options)
    res = strcmp(options.type,'rrt');
end

% R
function res = aux_getCondfunOptions_R(sys,func,params,options)
    res = any(strcmp(options.type,{'rrt','constrained'}));
end

% enclose
function res = aux_getCondfunOptions_enclose(sys,func,params,options)
    res = any(ismember(getMembers('guardIntersect4enclose'),options.guardIntersect));
end

% guardOrder
function res = aux_getCondfunOptions_guardOrder(sys,func,params,options)
    res = any(ismember(getMembers('guardIntersect4guardOrder'),options.guardIntersect));
end

% lagrangeRem.zooMethods
function res = aux_getCondfunOptions_lagrangeRem_zooMethods(sys,func,params,options)
    res = strcmp(options.lagrangeRem.method,'zoo');
end

% lagrangeRem.optMethod
function res = aux_getCondfunOptions_lagrangeRem_optMethod(sys,func,params,options)
    res = strcmp(options.lagrangeRem.method,'taylorModel');
end

% lagrangeRem.maxOrder
function res = aux_getCondfunOptions_lagrangeRem_maxOrder(sys,func,params,options)
    res = any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}));
end

% lagrangeRem.tolerance
function res = aux_getCondfunOptions_lagrangeRem_tolerance(sys,func,params,options)
    res = any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}));
end

% lagrangeRem.eps
function res = aux_getCondfunOptions_lagrangeRem_eps(sys,func,params,options)
    res = any(strcmp(options.lagrangeRem.method,{'taylorModel','zoo'}));
end

% intermediateOrder
function res = aux_getCondfunOptions_intermediateOrder(sys,func,params,options)
    if isa(sys,'nonlinParamSys') || strcmp(func,'reachInnerParallelotope')
        res = options.tensorOrder > 2;
    elseif strcmp(func,'reachInnerScaling')
        res = true;
    else
        res = ~contains(options.alg,'adaptive') && ~(strcmp(options.alg,'lin') && options.tensorOrder >= 2);
    end
end

% polyZono.maxDepGenOrder
function res = aux_getCondfunOptions_polyZono_maxDepGenOrder(sys,func,params,options)
    if startsWith(func,'reachInner')
        res = true;
    else
        res = contains(options.alg,'poly');
    end
end

% polyZono.maxPolyZonoRatio
function res = aux_getCondfunOptions_polyZono_maxPolyZonoRatio(sys,func,params,options)
    if startsWith(func,'reachInner')
        res = true;
    else
        res = contains(options.alg,'poly');
    end
end

% polyZono.restructureTechnique
function res = aux_getCondfunOptions_polyZono_restructureTechnique(sys,func,params,options)
    if startsWith(func,'reachInner')
        res = true;
    else
        res = contains(options.alg,'poly');
    end
end

% inpChanges
function res = aux_getCondfunOptions_inpChanges(sys,func,params,options)
    res = isfield(params,'U');
end


% ------------------------------ END OF CODE ------------------------------
