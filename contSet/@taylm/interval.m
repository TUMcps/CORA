function int = interval(tay,varargin)
% interval - Calculate an interval that bounds the taylor model
%
% Syntax:
%    int = interval(obj)
%    int = interval(obj,option)
%
% Inputs:
%    obj - taylm object
%    option - method used for the computation of the bounding interval
%             'int': standard interval arithmetic (default)
%             'bnb': branch and bound method is used to find min/max
%             'bnbAdv': branch and bound with re-expansion of taylor models
%             'linQuad': optimization with Linear Dominated Bounder (LDB)
%                        and Quadratic Fast Bounder (QFB) 
%             'bernstein': conversion to a bernstein polynomial
%
% Outputs:
%    int - interval overapproximating the taylm (class interval)
%
% Example:
%    tx = taylm(interval(1,4),4,'x');
%    ty = taylm(interval(1,4),4,'y');
%    tay = sin(tx + ty);
%    int1 = interval(tay,'int')
%    int2 = interval(tay,'bnb')
%
% Other m-files required: interval
% Subfunctions: s_tayl2int, intPower, intMul, evalInt
% MAT-files required: none
%
% See also: taylm, optBnb, optBnbAdv, optLinQuad
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"
%   [2] M. Althoff et al. "Implementation of Taylor models in CORA 2018
%       (Tool Presentation)"
%   [3] K. Makino et al. "Verified Global Optimization with Taylor Model 
%       based Range Bounders"

% Authors:       Niklas Kochdumper
% Written:       04-April-2018
% Last update:   ---
% Last revision: 17-August-2023 (TL, aux_ and comments)

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    option = setDefaultValues({tay.opt_method},varargin);

    % check input arguments
    inputArgsCheck({{tay,'att','taylm'};
                    {option,'str',{'int','bnb','bnbAdv','linQuad','bernstein'}}});    

    % calculate the bounding interval
    if strcmp(option, 'int')
        int = arrayfun(@(a) aux_tayl2int(a), tay, 'UniformOutput', false);
    elseif strcmp(option, 'bernstein')
        int = arrayfun(@(a) optBernstein(a), tay, 'UniformOutput', false);
    elseif strcmp(option, 'bnb')
        int = arrayfun(@(a) optBnb(a), tay, 'UniformOutput', false);
    elseif strcmp(option, 'bnbAdv')
        int = arrayfun(@(a) optBnbAdv(a), tay, 'UniformOutput', false);
    elseif strcmp(option, 'linQuad')
        int = arrayfun(@(a) optLinQuad(a), tay, 'UniformOutput', false);
    else
        throw(CORAerror('CORA:wrongValue','second',...
            'be ''int'', ''bnb'', ''bnbAdv'' and ''linQuad'''));
    end
    A = [int{:}];
    int = reshape(A, size(int));
    
end


% Auxiliary functions -----------------------------------------------------

function int = aux_tayl2int(tay)
    % evaluate taylor factors
    int = tay.remainder;
    
    for i = 1:length(tay.coefficients)
        exp = tay.monomials(i,:);
        temp = 1;       
        for j = 1:length(exp)
           temp = aux_intMul(temp,aux_intPower(exp(j)));
        end
        int = int + aux_evalInt(temp) * tay.coefficients(i);
    end
end

function int = aux_intPower(exponent)
    % exponent to integer

    if exponent == 0
        int = 1;            % interval(1,1)           
    elseif mod(exponent,2) == 0
        int = 2;            % interval(0,1)
    else
        int = 3;            % interval(-1,1)
    end
end

function int = aux_intMul(factor1,factor2)
    % also see aux_intPower

    if factor1 == 1
        int = factor2; % 1 * fac2 = fac2
    elseif factor2 == 1
        int = factor1; % fac1 * 1 = fac1
    elseif factor1 == 2 && factor2 == 2
        int = 2;       % [0,1] * [0,1] = [0,1]
    elseif factor1 == 3 && factor2 == 3
        int = 3;       % [-1,1] * [-1,1] = [-1,1]
    else
        int = 3;       % [-1,1] * [0,1] = [-1,1], [0,1] * [-1,1] = [-1,1]
    end
end

function int = aux_evalInt(tay)
    % also see aux_intPower
        
    if tay == 1
        int = interval(1,1);
    elseif tay == 2
        int = interval(0,1);
    else
        int = interval(-1,1);
    end
end

% ------------------------------ END OF CODE ------------------------------
