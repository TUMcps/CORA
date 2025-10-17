function checks = checkDynParameterOptions(field,sys,func,params,options,checks)
% checkDynParameterOptions - checks dynamic parameter values
%
% Syntax:
%    checkDynParameterOptions(field,sys,params,options,checks)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    func - function
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%    checks - struct
%
% Outputs:
%    checks - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkDynParameter

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for checks in options 
switch field % TODO sort and categorize
    case 'verbose'
        checks = aux_getChecksOptions_verbose(checks,sys,func,params,options);
    case 'reductionTechnique'
        checks = aux_getChecksOptions_reductionTechnique(checks,sys,func,params,options);
    case 'saveOrder'
        checks = aux_getChecksOptions_saveOrder(checks,sys,func,params,options);
        % alg
    case 'algInner'
        checks = aux_getChecksOptions_algInner(checks,sys,func,params,options);
    case 'reachAlg'
        checks = aux_getChecksOptions_reachAlg(checks,sys,func,params,options);
    case 'armaxAlg'
        checks = aux_getChecksOptions_armaxAlg(checks,sys,func,params,options);
    case 'linAlg'
        checks = aux_getChecksOptions_linAlg(checks,sys,func,params,options);
    case 'idAlg'
        checks = aux_getChecksOptions_idAlg(checks,sys,func,params,options);
    case 'alg'
        checks = aux_getChecksOptions_alg(checks,sys,func,params,options);
        % taylor
    case 'taylorOrder'
        checks = aux_getChecksOptions_taylorOrder(checks,sys,func,params,options);
    case 'taylmOrder'
        checks = aux_getChecksOptions_taylmOrder(checks,sys,func,params,options);
    case 'taylorTerms'
        checks = aux_getChecksOptions_taylorTerms(checks,sys,func,params,options);
        % ---
    case 'error'
        checks = aux_getChecksOptions_error(checks,sys,func,params,options);
    case 'timeStep'
        checks = aux_getChecksOptions_timeStep(checks,sys,func,params,options);
    case 'zonotopeOrder'
        checks = aux_getChecksOptions_zonotopeOrder(checks,sys,func,params,options);
    case 'l'
        checks = aux_getChecksOptions_l(checks,sys,func,params,options);
    case 'compOutputSet'
        checks = aux_getChecksOptions_compOutputSet(checks,sys,func,params,options);
    case 'partition'
        checks = aux_getChecksOptions_partition(checks,sys,func,params,options);
        % krylov
    case 'krylovError'
        checks = aux_getChecksOptions_krylovError(checks,sys,func,params,options);
    case 'krylovStep'
        checks = aux_getChecksOptions_krylovStep(checks,sys,func,params,options);
        % reduction and order
    case 'reductionTechniqueUnderApprox'
        checks = aux_getChecksOptions_reductionTechniqueUnderApprox(checks,sys,func,params,options);
    case 'tensorOrder'
        checks = aux_getChecksOptions_tensorOrder(checks,sys,func,params,options);
    case 'errorOrder'
        checks = aux_getChecksOptions_errorOrder(checks,sys,func,params,options);
        % ---
    case 'points'
        checks = aux_getChecksOptions_points(checks,sys,func,params,options);
    case 'p_conf'
        checks = aux_getChecksOptions_p_conf(checks,sys,func,params,options);
    case 'solver'
        checks = aux_getChecksOptions_solver(checks,sys,func,params,options);
    case 'vertSamp'
        checks = aux_getChecksOptions_vertSamp(checks,sys,func,params,options);
    case 'R'
        checks = aux_getChecksOptions_R(checks,sys,func,params,options);
    case 'type'
        checks = aux_getChecksOptions_type(checks,sys,func,params,options);
    case 'nrConstInp'
        checks = aux_getChecksOptions_nrConstInp(checks,sys,func,params,options);
    case 'stretchFac'
        checks = aux_getChecksOptions_stretchFac(checks,sys,func,params,options);
    case 'fracInpVert'
        checks = aux_getChecksOptions_fracInpVert(checks,sys,func,params,options);
    case 'fracVert'
        checks = aux_getChecksOptions_fracVert(checks,sys,func,params,options);
        % guard
    case 'guardIntersect'
        checks = aux_getChecksOptions_guardIntersect(checks,sys,func,params,options);
    case 'enclose'
        checks = aux_getChecksOptions_enclose(checks,sys,func,params,options);
    case 'guardOrder'
        checks = aux_getChecksOptions_guardOrder(checks,sys,func,params,options);
    case 'intersectInvariant'
        checks = aux_getChecksOptions_intersectInvariant(checks,sys,func,params,options);
    case 'compTimePoint'
        checks = aux_getChecksOptions_compTimePoint(checks,sys,func,params,options);
    case 'intermediateTerms'
        checks = aux_getChecksOptions_intermediateTerms(checks,sys,func,params,options);
    case 'gamma'
        checks = aux_getChecksOptions_gamma(checks,sys,func,params,options);
    case 'reductionInterval'
        checks = aux_getChecksOptions_reductionInterval(checks,sys,func,params,options);
        % max error
    case 'maxError'
        checks = aux_getChecksOptions_maxError(checks,sys,func,params,options);
    case 'maxError_x'
        checks = aux_getChecksOptions_maxError_x(checks,sys,func,params,options);
    case 'maxError_y'
        checks = aux_getChecksOptions_maxError_y(checks,sys,func,params,options);
    case 'tensorOrderOutput'
        checks = aux_getChecksOptions_tensorOrderOutput(checks,sys,func,params,options);
        % lagrangeRem
    case 'lagrangeRem.zooMethods'
        checks = aux_getChecksOptions_lagrangeRem_zooMethods(checks,sys,func,params,options);
    case 'lagrangeRem.simplify'
        checks = aux_getChecksOptions_lagrangeRem_simplify(checks,sys,func,params,options);
    case 'lagrangeRem.method'
        checks = aux_getChecksOptions_lagrangeRem_method(checks,sys,func,params,options);
    case 'lagrangeRem.tensorParallel'
        checks = aux_getChecksOptions_lagrangeRem_tensorParallel(checks,sys,func,params,options);
    case 'lagrangeRem.optMethod'
        checks = aux_getChecksOptions_lagrangeRem_optMethod(checks,sys,func,params,options);
    case 'lagrangeRem.replacements'
        checks = aux_getChecksOptions_lagrangeRem_replacements(checks,sys,func,params,options);
    case 'lagrangeRem.maxOrder'
        checks = aux_getChecksOptions_lagrangeRem_maxOrder(checks,sys,func,params,options);
    case 'lagrangeRem.tolerance'
        checks = aux_getChecksOptions_lagrangeRem_tolerance(checks,sys,func,params,options);
    case 'lagrangeRem.eps'
        checks = aux_getChecksOptions_lagrangeRem_eps(checks,sys,func,params,options);
        % ---
    case 'approxDepOnly'
        checks = aux_getChecksOptions_idDepOnly(checks,sys,func,params,options);
    case 'linearizationPoint'
        checks = aux_getChecksOptions_linearizationPoint(checks,sys,func,params,options);
    case 'intermediateOrder'
        checks = aux_getChecksOptions_intermediateOrder(checks,sys,func,params,options);
        % polyZono
    case 'polyZono.maxDepGenOrder'
        checks = aux_getChecksOptions_polyZono_maxDepGenOrder(checks,sys,func,params,options);
    case 'polyZono.maxPolyZonoRatio'
        checks = aux_getChecksOptions_polyZono_maxPolyZonoRatio(checks,sys,func,params,options);
    case 'polyZono.restructureTechnique'
        checks = aux_getChecksOptions_polyZono_restructureTechnique(checks,sys,func,params,options);
        % ---
    case 'timeStepInner'
        checks = aux_getChecksOptions_timeStepInner(charecks,sys,func,params,options);
    case 'contractor'
        checks = aux_getChecksOptions_contractor(checks,sys,func,params,options);
    case 'iter'
        checks = aux_getChecksOptions_iter(checks,sys,func,params,options);
    case 'splits'
        checks = aux_getChecksOptions_splits(checks,sys,func,params,options);
    case 'scaleFac'
        checks = aux_getChecksOptions_scaleFac(checks,sys,func,params,options);
    case 'orderInner'
        checks = aux_getChecksOptions_orderInner(checks,sys,func,params,options);
    case 'inpChanges'
        checks = aux_getChecksOptions_inpChanges(checks,sys,func,params,options);
    case 'approxErr'
        % ---
        checks = aux_getChecksOptions_idErr(checks,sys,func,params,options);
    case 'prevErrScale'
        checks = aux_getChecksOptions_prevErrScale(checks,sys,func,params,options);
    case 'prevErr'
        checks = aux_getChecksOptions_prevErr(checks,sys,func,params,options);
    case 'updateInitFnc'
        checks = aux_getChecksOptions_updateInitFnc(checks,sys,func,params,options);
    case 'norm'
        checks = aux_getChecksOptions_norm(checks,sys,func,params,options);
    case 'timeStepDivider'
        checks = aux_getChecksOptions_timeStepDivider(checks,sys,func,params,options);
    case 'postProcessingOrder'
        checks = aux_getChecksOptions_postProcessingOrder(checks,sys,func,params,options);
    %for conform:
    case 'cs.cost'
        checks = aux_getChecksOptions_cs_cost(checks,sys,func,params,options);
    case 'cs.constraints'
        checks = aux_getChecksOptions_cs_constraints(checks,sys,func,params,options);
    case 'cs.task'
        checks = aux_getChecksOptions_cs_task(checks,sys,func,params,options);
    case 'cs.numRotations'
        checks = aux_getChecksOptions_cs_numRotations(checks,sys,func,params,options);
        % scaling factors alpha
    case 'cs.a_min'
        checks = aux_getChecksOptions_cs_a_min(checks,sys,func,params,options);
    case 'cs.a_max'
        checks = aux_getChecksOptions_cs_a_max(checks,sys,func,params,options);
        % parameters
    case 'cs.set_p'
        checks = aux_getChecksOptions_cs_set_p(checks,sys,func,params,options);
    case 'cs.p0'
        checks = aux_getChecksOptions_cs_p0(checks,sys,func,params,options);
    case 'cs.cp_lim'
        checks = aux_getChecksOptions_cs_cp_lim(checks,sys,func,params,options);
        % optimization
    case 'cs.P'
        checks = aux_getChecksOptions_cs_P(checks,sys,func,params,options);
    case 'cs.verbose'
        checks = aux_getChecksOptions_cs_verbose(checks,sys,func,params,options);
    case 'cs.w'
        checks = aux_getChecksOptions_cs_w(checks,sys,func,params,options);
    case 'cs.robustness'
        checks = aux_getChecksOptions_cs_robustness(checks,sys,func,params,options);
    case 'cs.updateDeriv'
        checks = aux_getChecksOptions_cs_updateDeriv(checks,sys,func,params,options);
    case 'cs.timeout'
        checks = aux_getChecksOptions_cs_timeout(checks,sys,func,params,options);
        % outlier detection
    case 'cs.outMethod'
        checks = aux_getChecksOptions_cs_outMethod(checks,sys,func,params,options);
    case 'cs.numOutlier'
        checks = aux_getChecksOptions_cs_numOutlier(checks,sys,func,params,options);
    case 'cs.idzOutlier'
        checks = aux_getChecksOptions_cs_idzOutlier(checks,sys,func,params,options);
    case 'cs.determineActive'
        checks = aux_getChecksOptions_cs_determineActive(checks,sys,func,params,options);
        % 
    case 'cs.recMethod'
        checks = aux_getChecksOptions_cs_recMethod(checks,sys,func,params,options);
    case 'cs.batchSize'
        checks = aux_getChecksOptions_cs_batchSize(checks,sys,func,params,options);
    case 'cs.numCluster'
        checks = aux_getChecksOptions_cs_numCluster(checks,sys,func,params,options);
    case 'cs.numPoints'
        checks = aux_getChecksOptions_cs_numPoints(checks,sys,func,params,options);
    case 'cs.forgetting'
        checks = aux_getChecksOptions_cs_forgetting(checks,sys,func,params,options);
    case 'cs.forgettingCost'
        checks = aux_getChecksOptions_cs_forgettingCost(checks,sys,func,params,options);
    % identification
    case 'id.testSuite_val'
        checks = aux_getChecksOptions_id_testSuite_val(checks,sys,func,params,options);
    case 'id.p'
        checks = aux_getChecksOptions_id_p(checks,sys,func,params,options);
    case 'id.verbose'
        checks = aux_getChecksOptions_id_verbose(checks,sys,func,params,options);
    case 'id.filename'
        checks = aux_getChecksOptions_id_filename(checks,sys,func,params,options);
    case 'id.save_res'
        checks = aux_getChecksOptions_id_save_res(checks,sys,func,params,options);
        % genetic programming
    case 'id.gp_parallel'
        checks = aux_getChecksOptions_id_gp_parallel(checks,sys,func,params,options);
    case 'id.gp_runs'
        checks = aux_getChecksOptions_id_gp_runs(checks,sys,func,params,options);
    case 'id.gp_num_gen'
        checks = aux_getChecksOptions_id_gp_num_gen(checks,sys,func,params,options);
    case 'id.gp_max_genes'
        checks = aux_getChecksOptions_id_gp_max_genes(checks,sys,func,params,options);
    case 'id.gp_max_depth'
        checks = aux_getChecksOptions_id_gp_max_depth(checks,sys,func,params,options);
    case 'id.gp_func_names'
        checks = aux_getChecksOptions_id_gp_func_names(checks,sys,func,params,options);
    case 'id.gp_pop_size'
        checks = aux_getChecksOptions_id_gp_pop_size(checks,sys,func,params,options);
        % conformant genetic programming
    case 'id.cgp_num_gen'
        checks = aux_getChecksOptions_id_cgp_num_gen(checks,sys,func,params,options);
    case 'id.cgp_n_m_conf'
        checks = aux_getChecksOptions_id_cgp_n_m_conf(checks,sys,func,params,options);
    case 'id.cgp_pop_size_base'
        checks = aux_getChecksOptions_id_cgp_pop_size_base(checks,sys,func,params,options);
    otherwise
        CORAwarning('CORA:contDynamics','Unknown options.%s', field); return;
end


end


% Auxiliary functions -----------------------------------------------------

% options.<field> ---------------------------------------------------------

% verbose
function checks = aux_getChecksOptions_verbose(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% reductionTechnique
function checks = aux_getChecksOptions_reductionTechnique(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    if isa(sys,'nonlinearSys')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('reductionTechnique4nlsys'),val)), 'memberreductionTechnique4nlsys');
    elseif isfield(options,'linAlg') && startsWith(options.linAlg,'backward_minmax')
        % reductionTechniqueUnderApprox is used in this case
        checks = aux_getChecksOptions_reductionTechniqueUnderApprox(checks,sys,func,params,options);
    else
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('reductionTechnique'),val)), 'memberreductionTechnique');
    end
end

% saveOrder
function checks = aux_getChecksOptions_saveOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% linAlg
function checks = aux_getChecksOptions_linAlg(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    if isa(sys,'linearSysDT')
        % no member check
    elseif isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('linAlg4HA'),val)), 'memberlinAlg4HA');
    elseif contains(func,'Backward')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('linAlg4backward'),val)), 'memberlinAlg4backward');
    else
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('linAlg'),val)), 'memberlinAlg');
    end
end

% alg
function checks = aux_getChecksOptions_alg(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    if strcmp(func,'observe')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('alg4observe'),val)), 'memberalg4observe');
    elseif isa(sys,'nonlinDASys')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('alg4DA'),val)), 'memberalg4DA');
    elseif isa(sys,'nonlinearSysDT')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('alg4DT'),val)), 'memberalg4DT');
    elseif isa(sys,'nonlinearParamSys')
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('alg4param'),val)), 'memberalg4param');
    else
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('alg'),val)), 'memberalg');
    end
end

% error
function checks = aux_getChecksOptions_error(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% timeStep
function checks = aux_getChecksOptions_timeStep(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton') || isa(sys,'parallelHybridAutomaton')
        checks(end+1) = add2checks(@(val)c_HA_timeStep(val,sys,options), '');
    else
        checks(end+1) = add2checks(@isscalar, 'isscalar');
        checks(end+1) = add2checks(@(val)val>0, 'gezero');
        checks(end+1) = add2checks(@(val)abs((params.tFinal-params.tStart)/val - round((params.tFinal-params.tStart)/val))<1e-9, 'change_intsteps');
        checks(end+1) = add2checks(@(val)c_inputTraj(val,sys,params,options), '');
    end
end

% taylorTerms
function checks = aux_getChecksOptions_taylorTerms(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% zonotopeOrder
function checks = aux_getChecksOptions_zonotopeOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% l
function checks = aux_getChecksOptions_l(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)eq(size(val,1),sys.nrOfOutputs), 'eqoutput');
end

% compOutputSet
function checks = aux_getChecksOptions_compOutputSet(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% partition
function checks = aux_getChecksOptions_partition(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)c_partition(val,sys,options), '');
end

% krylovError
function checks = aux_getChecksOptions_krylovError(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% krylovStep
function checks = aux_getChecksOptions_krylovStep(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% reductionTechniqueUnderApprox
function checks = aux_getChecksOptions_reductionTechniqueUnderApprox(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('reductionTechniqueUnderApprox'),val)), 'memberreductionTechniqueUnderApprox');
end

% tensorOrder
function checks = aux_getChecksOptions_tensorOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)any(val==[2,3]), '2or3');
end

% errorOrder
function checks = aux_getChecksOptions_errorOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% points
function checks = aux_getChecksOptions_points(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% p_conf
function checks = aux_getChecksOptions_p_conf(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0) && le(val,1), 'normalized');
end

% solver
function checks = aux_getChecksOptions_solver(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
end

% type
function checks = aux_getChecksOptions_type(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)ismember(getMembers('type'),val), 'membertype');
end

% nrConstInp
function checks = aux_getChecksOptions_nrConstInp(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
    checks(end+1) = add2checks(@(val)c_nrConstInp(val,sys,params,options), '');
end

% fracInpVert
function checks = aux_getChecksOptions_fracInpVert(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0) && le(val,1), 'normalized');
end

% fracVert
function checks = aux_getChecksOptions_fracVert(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0) && le(val,1), 'normalized');
end

% vertSamp
function checks = aux_getChecksOptions_vertSamp(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% stretchFac
function checks = aux_getChecksOptions_stretchFac(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% R
function checks = aux_getChecksOptions_R(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)isa(val,'reachSet'), 'isareachSet');
end

% guardIntersect
function checks = aux_getChecksOptions_guardIntersect(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('guardIntersect'),val)), 'memberguardIntersect');
end

% enclose
function checks = aux_getChecksOptions_enclose(checks,sys,func,params,options)
    checks(end+1) = add2checks(@iscell, 'iscell');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('enclose'),val)), 'memberenclose');
end

% guardOrder
function checks = aux_getChecksOptions_guardOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% intersectInvariant
function checks = aux_getChecksOptions_intersectInvariant(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% compTimePoint
function checks = aux_getChecksOptions_compTimePoint(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% intermediateTerms
function checks = aux_getChecksOptions_intermediateTerms(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
    checks(end+1) = add2checks(@(val)le(val,options.taylorTerms), 'letaylorTerms');
end

% gamma
function checks = aux_getChecksOptions_gamma(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% reductionInterval
function checks = aux_getChecksOptions_reductionInterval(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
    checks(end+1) = add2checks(@(val)isinf(val)||mod(val,1)==0, 'integerorInf');
end

% maxError
function checks = aux_getChecksOptions_maxError(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isvector, 'isvector');
    checks(end+1) = add2checks(@(val)all(ge(val,0)), 'vectorgezero');
    checks(end+1) = add2checks(@(val)length(val)==sys.nrOfDims, 'eqsysdim');
end

% maxError_x
function checks = aux_getChecksOptions_maxError_x(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isvector, 'isvector');
    checks(end+1) = add2checks(@(val)all(ge(val,0)), 'vectorgezero');
    checks(end+1) = add2checks(@(val)length(val)==sys.nrOfDims, 'eqsysdim');
end

% maxError_y
function checks = aux_getChecksOptions_maxError_y(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isvector, 'isvector');
    checks(end+1) = add2checks(@(val)all(ge(val,0)), 'vectorgezero');
    checks(end+1) = add2checks(@(val)length(val)==sys.nrOfConstraints, 'eqconstr');
end

% tensorOrderOutput
function checks = aux_getChecksOptions_tensorOrderOutput(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)any(val==[2,3]), '2or3');
end

% lagrangeRem.simplify
function checks = aux_getChecksOptions_lagrangeRem_simplify(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('lagrangeRem.simplify'),val)), 'memberlagrangeRem.simplify');
end

% lagrangeRem.method
function checks = aux_getChecksOptions_lagrangeRem_method(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('lagrangeRem.method'),val)), 'memberlagrangeRem.method');
end

% lagrangeRem.tensorParallel
function checks = aux_getChecksOptions_lagrangeRem_tensorParallel(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% lagrangeRem.replacements
function checks = aux_getChecksOptions_lagrangeRem_replacements(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)isa(val,'function_handle'), 'isafunction_handle');
end

% lagrangeRem.zooMethods
function checks = aux_getChecksOptions_lagrangeRem_zooMethods(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('lagrangeRem.zooMethods'),val)), 'memberlagrangeRem.zooMethods');
end

% lagrangeRem.optMethod
function checks = aux_getChecksOptions_lagrangeRem_optMethod(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('lagrangeRem.optMethod'),val)), 'memberlagrangeRem.optMethod');
end

% lagrangeRem.maxOrder
function checks = aux_getChecksOptions_lagrangeRem_maxOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% lagrangeRem.tolerance
function checks = aux_getChecksOptions_lagrangeRem_tolerance(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% lagrangeRem.eps
function checks = aux_getChecksOptions_lagrangeRem_eps(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% approxDepOnly
function checks = aux_getChecksOptions_idDepOnly(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% linearizationPoint
function checks = aux_getChecksOptions_linearizationPoint(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isvector, 'isvector');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% intermediateOrder
function checks = aux_getChecksOptions_intermediateOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% polyZono.maxDepGenOrder
function checks = aux_getChecksOptions_polyZono_maxDepGenOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% polyZono.maxPolyZonoRatio
function checks = aux_getChecksOptions_polyZono_maxPolyZonoRatio(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% polyZono.restructureTechnique
function checks = aux_getChecksOptions_polyZono_restructureTechnique(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('restructureTechnique'),val)), 'memberrestructureTechnique');
end

% algInner
function checks = aux_getChecksOptions_algInner(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('algInner'),val)), 'memberalgInner');
end

% taylorOrder
function checks = aux_getChecksOptions_taylorOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% taylmOrder
function checks = aux_getChecksOptions_taylmOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% timeStepInner
function checks = aux_getChecksOptions_timeStepInner(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)abs(rem(params.tFinal,val))<eps, 'change_intsteps');
end

% contractor
function checks = aux_getChecksOptions_contractor(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('contractor'),val)), 'membercontractor');
end

% iter
function checks = aux_getChecksOptions_iter(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% splits
function checks = aux_getChecksOptions_splits(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% scaleFac
function checks = aux_getChecksOptions_scaleFac(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)c_scaleFac(val,sys,options), '');
end

% orderInner
function checks = aux_getChecksOptions_orderInner(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% inpChanges
function checks = aux_getChecksOptions_inpChanges(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
    checks(end+1) = add2checks(@(val)abs(rem(params.tFinal/(val+1),options.timeStep))<eps, 'eqreachSteps');
end

% approxErr
function checks = aux_getChecksOptions_idErr(checks,sys,func,params,options)
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% prevErrScale
function checks = aux_getChecksOptions_prevErrScale(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val) val>=0 && val<=1, 'normalized');
end

% updateInitFnc
function checks = aux_getChecksOptions_updateInitFnc(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)isa(val,'function_handle'), 'isafunction_handle');
end

% reachAlg
function checks = aux_getChecksOptions_reachAlg(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('reachAlg'),val)), 'memberreachAlg');
end

% norm
function checks = aux_getChecksOptions_norm(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('norm'),val)), 'membernorm');
end

% timeStepDivider
function checks = aux_getChecksOptions_timeStepDivider(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)val>0, 'gtzero');
end

% postProcessingOrder
function checks = aux_getChecksOptions_postProcessingOrder(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% armaxAlg
function checks = aux_getChecksOptions_armaxAlg(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('armaxAlg'),val)), 'memberarmaxAlg');
end

% idAlg
function checks = aux_getChecksOptions_idAlg(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('idAlg'),val)), 'memberidAlg');
end

% cs.cost
function checks = aux_getChecksOptions_cs_cost(checks,sys,func,params,options)
    %checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('cs.cost'),val)), 'membercs.cost');
end

% cs.constraints
function checks = aux_getChecksOptions_cs_constraints(checks,sys,func,params,options)
    %checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('cs.constraints'),val)), 'membercs.constraints');
end

% cs.task
function checks = aux_getChecksOptions_cs_task(checks,sys,func,params,options)
    %checks(end+1) = add2checks(@ischar, 'ischar');
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('cs.task'),val)), 'membercs.task');
end

% cs.numRotations
function checks = aux_getChecksOptions_cs_numRotations(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.set_p
function checks = aux_getChecksOptions_cs_set_p(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)isa(val,'function_handle'), 'isafunction_handle');
end

% cs.p0
function checks = aux_getChecksOptions_cs_p0(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% cs.cp_lim
function checks = aux_getChecksOptions_cs_cp_lim(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.a_min
function checks = aux_getChecksOptions_cs_a_min(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.a_max
function checks = aux_getChecksOptions_cs_a_max(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% cs.verbose
function checks = aux_getChecksOptions_cs_verbose(checks,sys,func,params,options)    
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% cs.robustness
function checks = aux_getChecksOptions_cs_robustness(checks,sys,func,params,options) 
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.updateDeriv
function checks = aux_getChecksOptions_cs_updateDeriv(checks,sys,func,params,options)    
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% cs.outMethod
function checks = aux_getChecksOptions_cs_outMethod(checks,sys,func,params,options) 
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('cs.outMethod'),val)), 'membercs.outMethod');
end

% cs.numOutlier
function checks = aux_getChecksOptions_cs_numOutlier(checks,sys,func,params,options) 
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.idzOutlier
function checks = aux_getChecksOptions_cs_idzOutlier(checks,sys,func,params,options) 
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val) all(val>= 1), 'vectorgeone');
    checks(end+1) = add2checks(@(val) all(mod(val,1)==0), 'integer');
end

% cs.determineActive
function checks = aux_getChecksOptions_cs_determineActive(checks,sys,func,params,options)    
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% cs.P
function checks = aux_getChecksOptions_cs_P(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% cs.w
function checks = aux_getChecksOptions_cs_w(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% cs.timeout
function checks = aux_getChecksOptions_cs_timeout(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% cs.recMethod
function checks = aux_getChecksOptions_cs_recMethod(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('cs.recMethod'),val)), 'membercs.recMethod');
end

% cs.batchSize
function checks = aux_getChecksOptions_cs_batchSize(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% cs.numPoints
function checks = aux_getChecksOptions_cs_numPoints(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% cs.numCluster
function checks = aux_getChecksOptions_cs_numCluster(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,1), 'geone');
end

% cs.forgetting
function checks = aux_getChecksOptions_cs_forgetting(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
    checks(end+1) = add2checks(@(val)le(val,1), 'leone');
end

% cs.forgettingCost
function checks = aux_getChecksOptions_cs_forgettingCost(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
    checks(end+1) = add2checks(@(val)le(val,1), 'leone');
end


% id.testSuite_val
function checks = aux_getChecksOptions_id_testSuite_val(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)isa(val,'trajectory'), 'istrajectory');
end

% id.p
function checks = aux_getChecksOptions_id_p(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% id.verbose
function checks = aux_getChecksOptions_id_verbose(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% id.filename
function checks = aux_getChecksOptions_id_filename(checks,sys,func,params,options)
    checks(end+1) = add2checks(@ischar, 'ischar');
end

% id.save_res
function checks = aux_getChecksOptions_id_save_res(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end

% id.nn_lr
function checks = aux_getChecksOptions_id_nn_lr(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.nn_lrSchedule
function checks = aux_getChecksOptions_id_nn_lrSchedule(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('id.nn_lrSchedule'),val)), 'memberid.nn_lrSchedule');
end

% id.nn_lrDropPeriod
function checks = aux_getChecksOptions_id_nn_lrDropPeriod(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.nn_lrDropFactor
function checks = aux_getChecksOptions_id_nn_lrDropFactor(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
    checks(end+1) = add2checks(@(val)le(val,1), 'leone');
end

% id.nn_act
function checks = aux_getChecksOptions_id_nn_act(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('id.nn_act'),val)), 'memberid.nn_act');
end

% id.nn_bs'
function checks = aux_getChecksOptions_id_nn_bs(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.nn_neurons
function checks = aux_getChecksOptions_id_nn_neurons(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
end

% id.nn_epochs
function checks = aux_getChecksOptions_id_nn_epochs(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_parallel
function checks = aux_getChecksOptions_id_gp_parallel(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@islogical, 'islogical');
end
 
% id.gp_runs
function checks = aux_getChecksOptions_id_gp_runs(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_num_gen
function checks = aux_getChecksOptions_id_gp_num_gen(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_max_genes
function checks = aux_getChecksOptions_id_gp_max_genes(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_max_depth
function checks = aux_getChecksOptions_id_gp_max_depth(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_func_names
function checks = aux_getChecksOptions_id_gp_func_names(checks,sys,func,params,options)
    add2checks(@iscell, 'iscell');
end

% id.cgp_num_gen
function checks = aux_getChecksOptions_id_cgp_num_gen(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.gp_pop_size
function checks = aux_getChecksOptions_id_gp_pop_size(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.cgp_n_m_conf
function checks = aux_getChecksOptions_id_cgp_n_m_conf(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% id.cgp_pop_size_base
function checks = aux_getChecksOptions_id_cgp_pop_size_base(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gtzero');
end

% ------------------------------ END OF CODE ------------------------------
