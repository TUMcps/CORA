function res = mpower(base,exponent)
% mpower - Overloaded '^' operator for intervals (power)
%
% Syntax:
%    res = mpower(base,exponent)
%
% Inputs:
%    base - interval object or numerical value
%    exponent - interval object or numerical value
%
% Outputs:
%    res - interval object
%
% Example:
%    base = interval([-2;1],[3;2]);
%    exponent = 2;
%    base^exponent
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Dmitry Grebenyuk
% Written:       25-June-2015
% Last update:   01-July-2015
%                10-February-2016
%                03-October-2017 (DG, A^0 = I is added)
% Last revision: 03-October-2017

% ------------------------------ BEGIN CODE -------------------------------

%scalar case
if isscalar(base) 
    if isscalar(exponent)
        res = base.^exponent;
        
        %see Appendix B of Introduction to Interval Analysis book
        %check whether exponent is an integer
        %if (abs(round(exponent) - exponent)) <= eps
        %    if exponent < 0
        %        res = 1/base^(-exponent);
        %    else
        %        if (base.inf > 0) || (mod(exponent,2))
        %            res.inf = base.inf^exponent;
        %            res.sup = base.sup^exponent;
        %        elseif (base.sup < 0) && ~(mod(exponent,2))
        %            res.inf = base.sup^exponent;
        %            res.sup = base.inf^exponent;
        %        elseif (base.inf <= 0 && base.sup >=0) && ~(mod(exponent,2))
        %            res.inf = 0;
        %            res.sup = max(abs(base.inf),abs(base.sup))^exponent;               
        %        end
        %    end
        %end
    end

%matrix case    
else 
    if isscalar(exponent) && exponent > 0
        %init
        res = base;
        for i = 2:exponent
            res = res*base;
        end
    elseif isscalar(exponent) && exponent == 0 && all(size(infimum(base)) == size(infimum(base)'))        
        res = eye(size(infimum(base),1));
    else
        throw(CORAerror('CORA:wrongValue','second','square matrix'));
    end
end

% ------------------------------ END OF CODE ------------------------------
