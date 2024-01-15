function res = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for intervals
%
% Syntax:
%    res = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - interval object
%    factor2 - interval object
%
% Outputs:
%    res - interval
%
% Example:
%    factor1 = interval(-2,0);
%    factor2 = interval(1,2);
%    factor1 * factor2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   25-June-2015
%                18-November-2015
%                01-February-2016 (DG, Fixed a matrix case)
%                27-February-2016 (DG, New matrix case)
%                21-July-2016 (MA, case scalar and numeric added)
%                22-July-2016 (MA, case that factor1 is numeric has been added)
%                26-July-2016 (multiplication with zonotope added)
%                05-August-2016 (simplified some cases; matrix case corrected)
%                23-June-2022 (VG, added support for empty matrices/intervals)
%                04-April-2023 (TL, vectorized matrix case)
%                10-October-2023 (TL, fix for scalar 0*[inf,-inf])
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%non-interval case
if isa(factor1,'interval') && ~isa(factor2,'interval')
    if isa(factor2,'zonotope') || isa(factor2,'conZonotope')
       res = intervalMultiplication(factor2,factor1);
       return;
    elseif isa(factor2,'polyZonotope') || isa(factor2,'zonoBundle') || ...
           isa(factor2,'conPolyZono')
       factor1 = intervalMatrix(center(factor1),rad(factor1));
       res = mtimes(factor1,factor2);
       return;
    end
end
    
%scalar case
if isscalar(factor1) && isscalar(factor2)
    
    %obtain possible values
    if isnumeric(factor1)
        res = factor2;
        if factor1 == 0
            % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
            possibleValues = 0;
        else
            possibleValues = [factor1*factor2.inf, factor1*factor2.sup];
        end
    elseif isnumeric(factor2)
        res = factor1;

        if factor2 == 0
            % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
            possibleValues = 0;
        else
            possibleValues = [factor1.inf*factor2, factor1.sup*factor2];
        end
    else
        res = factor1;
        possibleValues = [factor1.inf*factor2.inf, factor1.inf*factor2.sup, ...
            factor1.sup*factor2.inf, factor1.sup*factor2.sup];
    end
    
    %infimum
    res.inf = min(possibleValues);

    %supremum
    res.sup = max(possibleValues);

%mixed scalar/matric case: factor 1 is scalar
elseif isscalar(factor1)
    
    %obtain possible values
    if isnumeric(factor1) 
        if factor1 < 0
            %infimum and supremum
            res = interval( ...
                factor1*factor2.sup, ...
                factor1*factor2.inf ...
            );
        elseif factor1 > 0
            %infimum and supremum
            res = interval( ...
                factor1*factor2.inf, ...
                factor1*factor2.sup ...
            );
        else % factor1 == 0
            % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
            res = interval(zeros(size(factor2.inf)));
        end
    else
        res = factor1.*factor2;
    end

%mixed scalar/matrix case: factor 2 is scalar
elseif isscalar(factor2) 
    
    %obtain possible values
    if isnumeric(factor2)
        if factor2 < 0
            %infimum and supremum
            res = interval( ...
                factor2*factor1.sup, ...
                factor2*factor1.inf ...
            );
        elseif factor2 > 0
            %infimum and supremum
            res = interval( ...
                factor2*factor1.inf, ...
                factor2*factor1.sup ...
            );
        else % factor2 == 0
            % as 0*[-inf,inf] = {0*x|-inf<x<inf}={0}
            res = interval(zeros(size(factor1.inf)));
        end
    else
        res = factor1.*factor2;
    end

% matrix case
else
    if ~issparse(factor1) && ~issparse(factor2)
        % compute fast algorithm
        % [m, k] * [k, n] = [m, n]
        % -> [m, k, 1] .* [1, k, n] = [m, k, n]

        if isnumeric(factor2) 
            [m, ~] = size(factor1);
            extSize = [1, size(factor2)];
            factor2 = reshape(factor2, extSize);
            res = factor1 .* factor2;
        
            % [m,k,n] -> [m,n]
            res.inf = sum(res.inf, 2);
            res.sup = sum(res.sup, 2);
            res = reshape(res, m, []);
            
        else
            [m, ~] = size(factor1);
            extSize = [1, size(factor2)];
            factor2.inf = reshape(factor2.inf, extSize);
            factor2.sup = reshape(factor2.sup, extSize);
            res = factor1 .* factor2;
        
            % [m,k,n] -> [m,n]
            res.inf = sum(res.inf, 2);
            res.sup = sum(res.sup, 2);
            res = reshape(res, m, []);

        end
    else
        % sparse only supports 2d arrays, keeping old algorithm..

        % convert both to interval
        factor1 = interval(factor1);
        factor2 = interval(factor2);

        I1 = factor1.inf;
        S1 = factor1.sup;
    
        [m, ~] = size(I1);
        [~, n] = size(factor2.inf);
    
        % preallocate output bounds [m,n]
        resInf = zeros(m,n);
        resSup = zeros(m,n);
        
        A = interval(zeros(m,1));
        for i = 1:m
            % get i-th row [1, k], transpose
            A.inf = I1(i, :)';
            A.sup = S1(i, :)';

            % [k, 1] .* [k, n] = [k, n]
            B = A .* factor2;
            resInf(i, :) = sum(B.inf, 1);
            resSup(i, :) = sum(B.sup, 1);
        end
        res = interval(resInf, resSup);
    end
end


% ------------------------------ END OF CODE ------------------------------
