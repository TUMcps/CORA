function res = plus(factor1, factor2)
% plus - Overloaded '+' operator for a Taylor model
%
% Syntax:  
%    res = plus(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - a taylm objects
%    order  - the cut-off order of the Taylor series. The constat term is
%    the zero order.
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, minus, mtimes
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Dmitry Grebenyuk
% Written:      20-April-2016
%               21-July-2016 (DG) the polynomial part is changed to syms
%               18-July-2017 (DG) Multivariable polynomial pack is added
%               21-August-2017 (DG) Implementation fot matrices
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    if isscalar(factor1) && ~isscalar(factor2)
        res = arrayfun(@(b) s_plus(factor1, b), factor2, 'UniformOutput', 0);
    elseif ~isscalar(factor1) && isscalar(factor2)
        res = arrayfun(@(a) s_plus(a, factor2), factor1, 'UniformOutput', 0);  
    else
        res = arrayfun(@(a, b) s_plus(a, b), factor1, factor2, 'UniformOutput', 0);  
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end

%% --------------- Implementation for a scalar --------------

function res = s_plus(factor1, factor2)
    if isa(factor1, 'taylm') && isa(factor2, 'taylm')
        
        % find the common variables 
        [factor1, factor2] = rescale_dim(factor1, factor2);
        res = factor1;

        % Addition
        res.coefficients = [factor1.coefficients(:); factor2.coefficients(:)];
        res.monomials = [factor1.monomials; factor2.monomials];

        % Merge the properties of the two taylor models
        res = mergeProperties(res,factor1,factor2);
        
        % Reduce number of terms of the resulting Taylor model
        [res,rest] = compress(res);

        % Calculate the interval remainder
        res.remainder = factor1.remainder + factor2.remainder + rest;

        % if no polynomial part left, create an interval
        if isempty(res.coefficients) % subject to change. Maybe is wiser to keep the taylm
            res = res.remainder;
        end

    elseif isa(factor1,'taylm') && isa(factor2,'double')

        res = addConst(factor1,factor2);
    
    elseif isa(factor1,'double') && isa(factor2,'taylm')
        
        res = addConst(factor2,factor1);
        
    elseif isa(factor1,'taylm') && isa(factor2,'interval')
        
        res = factor1;
        res.remainder = res.remainder + factor2;
        
    elseif isa(factor1,'interval') && isa(factor2,'taylm') 
        
        res = factor2;
        res.remainder = res.remainder + factor1;
        
    else
        
        error('Wrong input')
        
    end

end

%% Auxiliary functions

function res = addConst(obj,const)

    res = obj;
        
    if isempty(obj.monomials) || obj.monomials(1,1) ~= 0
       res.coefficients = [const; obj.coefficients];
       res.monomials = [zeros(1, length(names_of_var) + 1); obj.monomials];
    else
       res.coefficients(1) = obj.coefficients(1) + const; 
    end

end