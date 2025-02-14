function res = rdivide(numerator, denominator)
% rdivide - Overload './' operator for Taylor models
%
% Syntax:
%    res = rdivide(numerator,denominator)
%
% Inputs:
%    numerator - numerator (class taylm)
%    denominator - denominator (class taylm)
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: priv_inverse
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Authors:       Niklas Kochdumper
% Written:       14-June-2017
% Last update:   ---  
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if isa(numerator,'taylm') && isa(denominator,'taylm')
        
        %res = numerator .* priv_inverse(denominator);
        
        res = arrayfun(@(a, b) a .* priv_inverse(b), numerator, denominator, 'UniformOutput', 0);
        
    
    elseif isa(numerator,'taylm') && isa(denominator,'double')
        
        %res = numerator .* (1/denominator);
        res = arrayfun(@(a) a .* (1/denominator), numerator, 'UniformOutput', 0);

        
    elseif isa(denominator,'taylm') && isa(numerator,'double')

        %res = numerator .* priv_inverse(denominator);
        res = arrayfun(@(b) numerator .* priv_inverse(b), denominator, 'UniformOutput', 0);

        
    end
    
    A = cat(1, res{:});
    res = reshape(A, size(res));
end

% ------------------------------ END OF CODE ------------------------------
