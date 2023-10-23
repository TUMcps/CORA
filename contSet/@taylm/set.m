function obj = set(obj,property,value)
% set - Set properties of a taylm object
%
% Syntax:
%    obj = set(obj,property,value)
%
% Inputs:
%    obj - taylm object
%    property - property that is changed ('max_order', 'opt_method', 'eps'
%               or 'tolerance')
%    value - new value of the property
%
% Outputs:
%    obj - taylm object
%
% Example:
%    tx = taylm(interval(1,4),4,'x');
%    ty = taylm(interval(1,4),4,'y');
%    t = sin(tx + ty);
%    t = set(t,'opt_method','int');
%    int1 = interval(t)
%    t = set(t,'opt_method','bnb');
%    int2 = interval(t)
%
% Other m-files required: interval
% Subfunctions: intPower, intMul, evalInt
% MAT-files required: none
%
% See also: taylm

% Authors:       Niklas Kochdumper
% Written:       06-April-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    obj = arrayfun(@(a) aux_s_set(a,property,value), obj, 'UniformOutput', 0);
    A = cat(1, obj{:});
    obj = reshape(A, size(obj));
    
end


% Auxiliary functions -----------------------------------------------------

function obj = aux_s_set(obj,property,value)

    if strcmp(property,'max_order')
        obj.max_order = value;
    elseif strcmp(property,'opt_method')
        obj.opt_method = value;
        if ~ismember(value,{'int','bnb','bnbAdv','linQuad'})
            throw(CORAerror('CORA:wrongValue','third',"'int', 'bnb', 'bnbAdv' or 'linQuad'"));
        end
    elseif strcmp(property,'eps')
        obj.eps = value;
    elseif strcmp(property,'tolerane')
        obj.tolerance = value;
    else
        throw(CORAerror('CORA:wrongValue','second',"'max_order', 'eps' or 'tolerance'"));
    end
end

% ------------------------------ END OF CODE ------------------------------
