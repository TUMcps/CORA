function res = times(factor1,factor2)
% times - Overloaded '.*' operator for intervals
%
% Syntax:
%    res = times(factor1,factor2)
%
% Inputs:
%    factor1 - interval object
%    factor2 - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    factor1 = interval([-2;1],[3;2]);
%    factor2 = interval([-1;-2],[1;2]);
%    factor1 .* factor2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   13-January-2016 (DG)
%                04-April-2023 (TL, minor optimizations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% an interval * a number
if isa(factor1, 'interval') && ~isa(factor2, 'interval')
    if isscalar(factor2)
        % use mtimes (considers scalar case explicitly)
        res = factor1 * factor2;
        return
    end

    resInf = factor1.inf .* factor2;
    resSup = factor1.sup .* factor2;

    res = factor1;
    res.inf = min(resInf, resSup);
    res.sup = max(resInf, resSup);
    
% a number * an interval
elseif ~isa(factor1, 'interval') && isa(factor2, 'interval')
    if isscalar(factor1)
        % use mtimes (considers scalar case explicitly)
        res = factor1 * factor2;
        return
    end
    
    resInf = factor2.inf .* factor1;
    resSup = factor2.sup .* factor1;
    
    res = factor2;
    res.inf = min(resInf, resSup);
    res.sup = max(resInf, resSup);
    
% an interval * an interval
else

    % possible combinations
    res1 = factor1.inf .* factor2.inf;
    res2 = factor1.inf .* factor2.sup;
    res3 = factor1.sup .* factor2.inf;
    res4 = factor1.sup .* factor2.sup;

    % to find min and max
    res = factor1;
    res.inf = min(res1,min(res2, min(res3,res4)));
    res.sup = max(res1,max(res2, max(res3,res4)));    
end

% ------------------------------ END OF CODE ------------------------------
