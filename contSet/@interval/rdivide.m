function res = rdivide(numerator, denominator)
% rdivide - Overloads the ./ operator that provides elementwise division
%    of two matrices
%
% For an interval / a number
% [-Inf, +Inf]                                      if denominator = 0;
% [min(inf / n, sup / n), max(inf / n, sup / n)     if denominator ~= 0.
%
% For a number / an interval
% [min(n / inf, n / sup), max(n / inf, n / sup)     if inf > 0 or sup < 0;
% [ n / sup, +Inf]                                  if inf = 0;
% [ -Inf, n / inf]                                  if sup = 0;
% [NaN, NaN]                                        if inf = 0 and sup = 0;
% [-Inf, +Inf]                                      if inf < 0 and sup > 0.
%
% Syntax:
%    res = rdivide(numerator, denominator)
%
% Inputs:
%    numerator - interval object
%    denominator - interval object
%
% Outputs:
%    res - interval object after elementwise division
%
% Example: 
%    I = interval([-4;2],[1;3]);
%    divisor = [3;2];
%    I = I./divisor
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Dmitry Grebenyuk
% Written:       07-February-2016
% Last update:   13-March-2016 (speed improvement)
%                05-May-2020 (MW, added error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% an interval / a (matrix of) scalar(s)
if isa(numerator, 'interval') && ~isa(denominator, 'interval')
    
    % a matrix of scalars or a scalar
    if isequal(size(denominator),size(numerator)) || isequal(size(denominator),[1, 1])

        % needs to preserve the shape of an input
        res = numerator;
        res.inf = min( numerator.inf ./ denominator, numerator.sup ./ denominator );
        res.sup = max( numerator.inf ./ denominator, numerator.sup ./ denominator );       

        if isequal(size(denominator), size(numerator))
            ind1 = find(denominator == 0);
            res.inf(ind1) = -Inf;
            res.sup(ind1) = +Inf;
        elseif denominator == 0
            res.inf(:) = -Inf;
            res.sup(:) = +Inf;
        end

    else
        throw(CORAerror('CORA:specialError','The input size is wrong.'));
    end

% a (matrix of) scalar(s) / an interval
elseif ~isa(numerator, 'interval') && isa(denominator, 'interval')
    
    if size(numerator) == 1

        % needs to preserve the shape of an input
        res = denominator;
        res.inf = min(numerator ./ denominator.sup, numerator ./ denominator.inf);
        res.sup = max(numerator ./ denominator.sup, numerator ./ denominator.inf);

        ind1 = denominator.inf == 0 & denominator.sup == 0;
        res.inf(ind1) = NaN;
        res.sup(ind1) = NaN;

        ind1 = denominator.sup == 0;
        res.inf(ind1) = -Inf;
        res.sup(ind1) = numerator ./ denominator.inf(ind1);

        ind1 = denominator.inf < 0 & denominator.sup > 0;
        res.inf(ind1) = -Inf;
        res.sup(ind1) = +Inf;   

    elseif size(denominator) == size(numerator) 

        % needs to preserve the shape of an input
        res = denominator;
        res.inf = min(numerator ./ denominator.sup, numerator ./ denominator.inf);
        res.sup = max(numerator ./ denominator.sup, numerator ./ denominator.inf);

        ind1 = find( denominator.inf == 0 & denominator.sup == 0);
        res.inf(ind1) = [];
        res.sup(ind1) = [];

        ind1 = find( denominator.sup == 0);
        res.inf(ind1) = -Inf;
        res.sup(ind1) = numerator(ind1) ./ denominator.inf(ind1);

        ind1 = find( denominator.inf < 0 & denominator.sup > 0);
        res.inf(ind1) = -Inf;
        res.sup(ind1) = +Inf;

    else
        throw(CORAerror('CORA:specialError','The input size is wrong.'));
    end
    
% an interval / an interval (x / y)
elseif isa(numerator, 'interval') && isa(denominator, 'interval')
    
    if isequal(size(numerator), size(denominator)) || ...
            isequal(size(numerator), [1, 1]) || ...
            isequal(size(denominator), [1, 1])
        y = 1 ./ denominator;
        res = numerator .* y;
    else
        throw(CORAerror('CORA:specialError','The input size is wrong.'));
    end
    
end

% return error if NaN occures
if any(any(isnan(res.inf))) || any(any(isnan(res.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','inf ~= 0 && sup ~= 0'));
end

% ------------------------------ END OF CODE ------------------------------
