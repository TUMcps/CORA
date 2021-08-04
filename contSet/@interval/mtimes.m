function res = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for intervals
%
% Syntax:  
%    res = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - interval (for computational efficiency, no single value
%              considered; does not require type checking)
%    factor2 - interval (for computational efficiency, no single value
%              considered; does not require type checking)
%
% Outputs:
%    res - interval
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      19-June-2015
% Last update:  25-June-2015
%               18-November-2015
%               01-February-2016 (DG, Fixed a matrix case)
%               27-February-2016(DG, New matrix case)
%               21-July-2016 (MA, case that factor1 is a scalar interval
%                                 and factor 2 is numeric added)
%               22-July-2016 (MA, case that factor1 is numeric has been added)
%               26-July-2016 (multiplication with zonotope added)
%               05-August-2016 (simplified some cases; matrix case corrected)
% Last revision:---

%------------- BEGIN CODE --------------

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
        possibleValues = [factor1*factor2.inf, factor1*factor2.sup];
    elseif isnumeric(factor2)
        res = factor1;
        possibleValues = [factor1.inf*factor2, factor1.sup*factor2];
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
        res = factor2;
        if factor1<0
            %infimum and supremum
            res.inf = factor1*factor2.sup;
            res.sup = factor1*factor2.inf;
        else
            %infimum and supremum
            res.inf = factor1*factor2.inf;
            res.sup = factor1*factor2.sup;
        end
    else
        res = factor1.*factor2;
    end

%mixed scalar/matrix case: factor 2 is scalar
elseif isscalar(factor2) 
    
    %obtain possible values
    if isnumeric(factor2)
        res = factor1;
        if factor2<0
            %infimum and supremum
            res.inf = factor2*factor1.sup;
            res.sup = factor2*factor1.inf;
        else
            %infimum and supremum
            res.inf = factor2*factor1.inf;
            res.sup = factor2*factor1.sup;
        end
    else
        res = factor1.*factor2;
    end

% matrix case [int] * [numb]

elseif isa(factor1, 'interval') && ~isa(factor2, 'interval')
    I1 = factor1.inf;
    S1 = factor1.sup;
 
    [m, n] = size(I1);
    [m1, n1] = size(factor2);
    A = interval();
 
    for i = 1:m
        A.inf = repmat(I1(i, :),n1, 1)';
        A.sup = repmat(S1(i, :),n1, 1)';
        B = A .* factor2;
        Binf(i, :) = sum(B.inf, 1);
        Bsup(i, :) = sum(B.sup, 1);
    end
    B.inf = Binf;
    B.sup = Bsup;
 
    res = B;
     
% matrix case [numb] * [int]

elseif ~isa(factor1, 'interval') && isa(factor2, 'interval')
    I1 = factor2.inf;
    S1 = factor2.sup;

    [m1, n1] = size(I1);
    [m, n] = size(factor1);
    A = interval();

    for i = 1:m
%         A.inf = repmat(I1(i, :),n1, 1);
%         A.sup = repmat(S1(i, :),n1, 1);
        factor1_1 = repmat(factor1(i, :), n1, 1)';

        B = factor1_1 .* factor2;
        Binf(i, :) = sum(B.inf, 1);
        Bsup(i, :) = sum(B.sup, 1);
    end
    B.inf = Binf;
    B.sup = Bsup;

    res = B;

% matrix case [int] * [int]
else
    
    I1 = factor1.inf;
    S1 = factor1.sup;

    [m, n] = size(I1);
    [m1, n1] = size(factor2.inf);
    A = interval();

    for i = 1:m
        A.inf = repmat(I1(i, :),n1, 1)';
        A.sup = repmat(S1(i, :),n1, 1)';
        B = A .* factor2;
        Binf(i, :) = sum(B.inf, 1);
        Bsup(i, :) = sum(B.sup, 1);
    end
    B.inf = Binf;
    B.sup = Bsup;

    res = B;

end


%------------- END OF CODE --------------