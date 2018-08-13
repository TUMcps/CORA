function int = interval(obj,varargin)
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
%
% Outputs:
%    int - interval overapproximating the taylm (class interval)
%
% Example:
%    tx = taylm(interval(1,4),4,'x');
%    ty = taylm(interval(1,4),4,'y');
%    t = sin(tx + ty);
%    int1 = interval(t,'int')
%    int2 = interval(t,'bnb')
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

% Author:       Niklas Kochdumper
% Written:      04-April-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE -------------

    % parse input
    option = obj.opt_method;
    if nargin == 2
       option = varargin{1};
    end

    % calculate the bounding interval
    if strcmp(option, 'int')
        int = arrayfun(@(a) s_tayl2int(a), obj, 'UniformOutput', 0);
    elseif strcmp(option, 'bnb')
        int = arrayfun(@(a) optBnb(a), obj, 'UniformOutput', 0);
    elseif strcmp(option, 'bnbAdv')
        int = arrayfun(@(a) optBnbAdv(a), obj, 'UniformOutput', 0);
    elseif strcmp(option, 'linQuad')
        int = arrayfun(@(a) optLinQuad(a), obj, 'UniformOutput', 0); 
    else
        error('Only values ''int'', ''bnb'', ''bnbAdv'' and ''linQuad'' are valid for input argument ''option''!');
    end
    A = [int{:}];
    int = reshape(A, size(int));
    
end


%% Auxiliary functions

function int = s_tayl2int(obj)
    % evaluate taylor factors
    int = obj.remainder;
    
    for i = 1:length(obj.coefficients)
        exp = obj.monomials(i,:);
        temp = 1;       
        for j = 1:length(exp)
           temp = intMul(temp,intPower(exp(j)));
        end
        int = int + evalInt(temp) * obj.coefficients(i);
    end

end

function int = intPower(exponent)

    if exponent == 0
        int = 1;            % interval(1,1)           
    elseif mod(exponent,2) == 0
        int = 2;            % interval(0,1)
    else
        int = 3;            % interval(-1,1)
    end
    
end

function int = intMul(factor1,factor2)

    if factor1 == 1
        int = factor2;
    elseif factor2 == 1
        int = factor1;
    elseif factor2 == 2 && factor1 == 2
        int = 2;      % [0,1] * [0,1] = [0,1]
    elseif factor2 == 3 && factor2 == 3
        int = 3;      % [-1,1] * [-1,1] = [-1,1]
    else
        int = 3;      % [-1,1] * [0,1] = [-1,1]
    end
    
end

function int = evalInt(obj)
        
    if obj == 1
        int = interval(1,1);
    elseif obj == 2
        int = interval(0,1);
    else
        int = interval(-1,1);
    end

end

%------------ END OF CODE ------------ 