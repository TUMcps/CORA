function [objX,objU] = initRangeBoundingObjects(intX,intU,options)
% initRangeBoundingObjects - creates taylm- or zoo-objects for the state
%    and input variables
%
% Syntax:
%    [objX, objU] = initRangeBoundingObjects(intX, intU, options)
%
% Inputs:
%    intX - interval bounding the state variables
%    intU - interval bounding the input variables
%    options - struct containing algorithm settings
%
% Outputs:
%    objX - taylm-/zoo-object for the state variables
%    objU - taylm-/zoo-object for the input variables
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, zoo

% Authors:       Niklas Kochdumper
% Written:       02-August-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse options
    maxOrder = [];
    optMethod = [];
    eps = [];
    tolerance = [];

    if isfield(options.lagrangeRem,'maxOrder')
        maxOrder = options.lagrangeRem.maxOrder;
    end
    if isfield(options.lagrangeRem,'optMethod')
        optMethod = options.lagrangeRem.optMethod;
    end
    if isfield(options.lagrangeRem,'tolerance')
        eps = options.lagrangeRem.eps;
    end
    if isfield(options.lagrangeRem,'eps')
        tolerance = options.lagrangeRem.tolerance;
    end

    % generate taylor model or zoo objects
    if strcmp(options.lagrangeRem.method,'taylorModel')

        objX = taylm(intX,maxOrder,aux_idxVars('x',1:length(intX)), ...
                     optMethod,eps,tolerance);
        objU = taylm(intU,maxOrder,aux_idxVars('u',1:length(intU)), ...
                     optMethod,eps,tolerance);

    elseif strcmp(options.lagrangeRem.method,'zoo')

        objX = zoo(intX,options.lagrangeRem.zooMethods, ...
                   aux_idxVars('x',1:length(intX)),maxOrder,eps,tolerance);
        objU = zoo(intU,options.lagrangeRem.zooMethods, ...
                   aux_idxVars('u',1:length(intU)),maxOrder,eps,tolerance);

    else
        throw(CORAerror('CORA:wrongFieldValue','options.lagrangeRem.method',...
            {'taylorModel','zoo'}));
    end

end


% Auxiliary functions -----------------------------------------------------

function indexedVars = aux_idxVars(varName,idx)

    idxLength = length(idx);
    indexedVars = cell(idxLength,1);
    for i=1:idxLength
        indexedVars{i} = [varName,num2str(idx(i))];
    end

end

% ------------------------------ END OF CODE ------------------------------
