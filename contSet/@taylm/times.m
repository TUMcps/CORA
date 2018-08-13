function res = times(factor1,factor2)
% times - Overload '.*' operator for Taylor models
%
% Syntax:  
%    res = times(factor1,factor2)
%
% Inputs:
%    factor1 - first taylm object
%    factor2 - second taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval, interval
% Subfunctions: multiply
% MAT-files required: none
%
% See also: taylm, plus, minus
%
% References: 
%   [1] K. Makino et al. "Taylor Models and other validated functional 
%       inclusion methods"

% Author:       Niklas Kochdumper, Dmitry Grebenyuk
% Written:      14-June-2017
%               11-November-2017 (DG) extra cases are added
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------
    
    if isa(factor1,'taylm') && isa(factor2,'taylm') && all(size(factor1) == size(factor2))
        
        res = arrayfun(@(a, b) s_times_tt(a, b), factor1, factor2, 'UniformOutput', 0);
        
    elseif isa(factor1,'taylm') && isa(factor2,'double') && isscalar(factor2)
        
        res = arrayfun(@(a) s_times_td(a, factor2), factor1, 'UniformOutput', 0);
        
    elseif isa(factor1,'taylm') && isscalar(factor1) && isa(factor2,'double')
        
        res = arrayfun(@(b) s_times_td(factor1, b), factor2, 'UniformOutput', 0);
        
    elseif isa(factor1,'taylm') && isa(factor2,'double') && all(size(factor1) == size(factor2))
        
        res = arrayfun(@(a, b) s_times_td(a, b), factor1, factor2, 'UniformOutput', 0);
        
    elseif isa(factor1,'double') && isscalar(factor1) && isa(factor2,'taylm') 

        res = arrayfun(@(b) s_times_dt(factor1, b), factor2, 'UniformOutput', 0);
        
    elseif isa(factor1,'double') && isa(factor2,'taylm') && isscalar(factor2)

        res = arrayfun(@(a) s_times_dt(a, factor2), factor1, 'UniformOutput', 0);
        
    elseif isa(factor1,'double') && isa(factor2,'taylm') && all(size(factor1) == size(factor2))

        res = arrayfun(@(a, b) s_times_dt(a, b), factor1, factor2, 'UniformOutput', 0);
        
    %elseif isa(factor1,'taylm') && isa(factor2,'interval') && all(size(factor1) == size(factor2))
        
    %    res = arrayfun(@(a, b) s_times_ti(a, b), factor1, factor2, 'UniformOutput', 0);
        
    %elseif isa(factor1,'taylm') && isscalar(factor1) && isa(factor2,'interval') 
        
    %    res = arrayfun(@(b) s_times_ti(factor1, b), factor2,'UniformOutput', 0);
        
    %elseif isa(factor1,'taylm') && isa(factor2,'interval') && isscalar(factor2)
        
    %    res = arrayfun(@(a) s_times_ti(a, factor2), factor1,'UniformOutput', 0);
        
    %elseif isa(factor1,'interval') && isa(factor2,'taylm') && all(size(factor1) == size(factor2))
        
    %    res = arrayfun(@(a, b) s_times_it(a, b), factor1, factor2, 'UniformOutput', 0);
        
    %elseif isa(factor1,'interval') && isscalar(factor1) && isa(factor2,'taylm')
        
    %    res = arrayfun(@(b) s_times_it(factor1, b), factor2, 'UniformOutput', 0);
        
    %elseif isa(factor1,'interval') && isa(factor2,'taylm') && isscalar(factor2)
        
    %    res = arrayfun(@(a) s_times_it(a, factor2), factor1, 'UniformOutput', 0);
        
    else
        
        error('Wrong input')
        
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end

%% --------------- Implementation for a scalar --------------
function res = s_times_tt(factor1, factor2)
            
    [factor1, factor2] = rescale_dim(factor1, factor2);
    res = factor1;

    % Multiplication
    [res.coefficients, res.monomials] = multiply ( ...
                        factor1.coefficients, factor1.monomials, ...
                        factor2.coefficients, factor2.monomials);

    % Merge the properties of the two taylor models
    res = mergeProperties(res,factor1,factor2);
                    
    % Reduce number of terms of resulting Taylor Model
    [res,rest] = compress(res);

    % Calculate remainder
    remainder1 = factor1.remainder;
    remainder2 = factor2.remainder;
    factor1.remainder = interval(0,0);
    factor2.remainder = interval(0,0);

    rem1 = rest;
    rem2 = interval(factor1) * remainder2;
    rem3 = interval(factor2) * remainder1;

    rem4 = remainder1 * remainder2;

    res.remainder = rem1 + rem2 + rem3 + rem4; 
    
end

function res = s_times_td(factor1, factor2)

    res = factor1;
    res.coefficients = res.coefficients * factor2;
    res.remainder = res.remainder * factor2;
    
end

function res = s_times_dt(factor1, factor2)

    res = factor2;
    res.coefficients = res.coefficients * factor1;
    res.remainder = res.remainder * factor1;
    
end

function res = s_times_ti(factor1, factor2)

    tolerance = 1e-8;

    % interpret intervals like [2,2] as a number
    if abs( supremum(factor2) - infimum(factor2) ) <= tolerance
        res = factor1 .* mid(factor2);

    else

        res = factor1 .* taylm(factor2, factor1.max_order,...
                            {strcat('inpn', num2str(round(rand(1)*100)) ) } );
                        % inputname(2) returns factor2 but no name of
                        % a variable, hence quick fix. To fix in
                        % future.
    end
    
end

function res = s_times_it(factor1, factor2)

    res = factor2 .* factor1;

end
%% Auxiliary functions

function [c,e] = multiply(c1,e1,c2,e2)

  [len, width] = size(e1);
  len =  len * length(c2);
  
  e = zeros ( len, width );
  c = zeros( len, 1 );

  len = 0;
  for j = 1 : length(c2)
    for i = 1 : length(c1)
      len = len + 1;
      c(len) = c1(i) * c2(j);
      e(len,:) = e1(i,:) + e2(j,:);
    end
  end

end

%------------ END OF CODE ------------
    