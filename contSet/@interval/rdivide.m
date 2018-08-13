function res = rdivide(numerator, denominator)
% rdivide - Overloads the ./ operator that provides elementwise division of 
% two matrices
%
% Syntax:  
%    res = rdivide( numerator, denominator )
%
% Inputs:
%    numerator, denominator - interval objects
%
% Outputs:
%    res - interval object after elementwise division
%
% Example: 
%    IH=intervalhull(rand(6,2));
%    divisor=rand(6,1);
%    IH=IH./divisor
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author:       Dmitry Grebenyuk
% Written:      07-February-2016
% Last update:  13-March-2016       Speed improvement
% Last revision:---

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

%------------- BEGIN CODE --------------

% an interval / a (maxtrix of) scalar(s)
if isa(numerator, 'interval') == 1 && isa(denominator, 'interval') == 0
    
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
               error('The input size is wrong.')
           end
            

% a (matrix of) scalar(s) / an interval
elseif isa(numerator, 'interval') == 0 && isa(denominator, 'interval') == 1
    
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
               error('The input size is wrong.')
            end
    
% an interval / an interval (x / y)
elseif isa(numerator, 'interval') == 1 && isa(denominator, 'interval') == 1
    
    if isequal(size(numerator), size(denominator)) || isequal(size(numerator), [1, 1]) == 1 || isequal(size(denominator), [1, 1]) == 1
        y = 1 ./ denominator;
        res = numerator .* y;
    else
        error('The input size is wrong.')
    end
    
end

%------------- END OF CODE --------------