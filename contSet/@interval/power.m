function res = power(base,exponent)
% power - Overloaded '.^' operator for intervals (power)
%
% For an interval .^ an integer number
% [min(base.^exp), max(base.^exp)   if the exp is an integer and odd;
% [0, max(base.^exp)                if the exp is an integer and even.
%
% For an interval .^ a real number 
% [min(base.^exp), max(base.^exp)   if base.inf >= 0;
% [NaN, NaN]                        if otherwise.
%
% For a number .^ an interval
% [min(base.^exp), max(base.^exp)   if base.inf >= 0;
% [NaN, NaN]                        if otherwise.
%
% For an interval .^ an interval
% [min(base.^exp), max(base.^exp)   if base.inf >= 0;
% [NaN, NaN]                        if otherwise.
%
% Syntax:  
%    res = power(base,exponent)
%
% Inputs:
%    base - interval object or numerical value
%    exponent - interval object or numerical value
%
% Outputs:
%    res - interval object
%
% Example:
%    base = interval([-3;2],[5;4]);
%    exponent = [2;3];
%    base .^ exponent
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Dmitry Grebenyuk
% Written:      10-February-2016
% Last update:  15-March-2016 (DG, speed up)
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% an interval .^ a number
if isnumeric(exponent)
    
    res = base;

    if isscalar(exponent)
       
        if (abs(round(exponent) - exponent)) <= eps
            
            % positive scalar integer exponent
            if exponent >= 0

                temp1 = base.inf .^ exponent;
                temp2 = base.sup .^ exponent;
                res.inf = min(temp1, temp2);
                res.sup = max(temp1, temp2);   

                % special behavior for even exponents
                if ~(mod(exponent,2)) && exponent ~= 0
                    ind1 =  base.inf < 0 & base.sup > 0 & ~(mod(exponent,2)) & exponent ~= 0;
                    res.inf(ind1) = 0;
                end
                
            % negative scalar integer exponent
            else
                res = (1 ./ base) .^ -exponent;
            end 
            
        else
            
            % positive scalar real valued exponent
            if exponent >= 0
               
                if isscalar(base)
                    if base.inf >= 0
                        res.inf = base.inf .^ exponent;
                        res.sup = base.sup .^ exponent;
                    else
                        res.inf = NaN;
                        res.sup = NaN;
                    end
                    
                else
                    
                    res.inf = base.inf .^ exponent;
                    res.sup = base.sup .^ exponent;
                    
                    ind1 = base.inf < 0;
                    res.inf(ind1) = NaN;
                    res.sup(ind1) = NaN;                 
                end 
                
            % negative scalar real valued exponent
            else                
                res = (1 ./ base) .^ -exponent;                
            end          
        end
 
    else
        
        % matrix exponent with scalar base
        if isscalar(base)
            
            res.inf = zeros(size(exponent));
            res.sup = zeros(size(exponent));
            
            % all integer exponent matrix
            if all(all((abs(round(exponent) - exponent)) <= eps))

                ind1 = exponent < 0;
                if any(any(ind1))
                    oneover = 1 ./ base;
                    temp1 = oneover.inf .^ -exponent(ind1);
                    temp2 = oneover.sup .^ -exponent(ind1);
                    res.inf(ind1) = min(temp1, temp2);
                    res.sup(ind1) = max(temp1, temp2);
                end

                ind1 = exponent >= 0;
                temp1 = base.inf .^ exponent(ind1);
                temp2 = base.sup .^ exponent(ind1);
                res.inf(ind1) = min(temp1, temp2);
                res.sup(ind1) = max(temp1, temp2);

                if base.inf < 0 && base.sup > 0
                    ind1 = ~(mod(exponent,2)) & exponent ~= 0;
                    res.inf(ind1) = 0; 
                end

            % mixed real and integer exponent matrix    
            else 
                
                ind1 = exponent < 0;
                if any(any(ind1))
                    oneover = 1 ./ base;
                    temp1 = oneover.inf .^ -exponent(ind1);
                    temp2 = oneover.sup .^ -exponent(ind1);
                    res.inf(ind1) = min(temp1, temp2);
                    res.sup(ind1) = max(temp1, temp2);
                end

                ind1 = exponent >= 0;
                temp1 = base.inf .^ exponent(ind1);
                temp2 = base.sup .^ exponent(ind1);
                res.inf(ind1) = min(temp1, temp2);
                res.sup(ind1) = max(temp1, temp2);
                
                if base.inf < 0 && base.sup > 0
                    ind1 = ~(mod(exponent,2)) & exponent ~= 0;
                    res.inf(ind1) = 0; 
                end
                
                if base.inf < 0
                    ind1 = abs(round(exponent)-exponent) > eps;
                    res.inf(ind1) = NaN;
                    res.sup(ind1) = NaN;
                end
            end
            
            
        % matrix exponent with matrix base
        else
            
            if all(size(base) == size(exponent))
                
                % all integer exponent matrix
                if all(all((abs(round(exponent) - exponent)) <= eps))
                
                    ind1 = exponent < 0;
                    if any(any(ind1))
                        oneover = 1 ./ base;
                        temp1 = oneover.inf(ind1) .^ -exponent(ind1);
                        temp2 = oneover.sup(ind1) .^ -exponent(ind1);
                        res.inf(ind1) = min(temp1, temp2);
                        res.sup(ind1) = max(temp1, temp2);
                    end

                    ind1 = exponent >= 0;
                    temp1 = base.inf(ind1) .^ exponent(ind1);
                    temp2 = base.sup(ind1) .^ exponent(ind1);
                    res.inf(ind1) = min(temp1, temp2);
                    res.sup(ind1) = max(temp1, temp2);

                    ind1 = base.inf < 0 & base.sup > 0 & ~(mod(exponent,2)) & exponent ~= 0;
                    res.inf(ind1) = 0; 
                    
                % mixed real and integer exponent matrix    
                else
                    
                    ind1 = exponent < 0;
                    if any(any(ind1))
                        oneover = 1 ./ base;
                        temp1 = oneover.inf(ind1) .^ -exponent(ind1);
                        temp2 = oneover.sup(ind1) .^ -exponent(ind1);
                        res.inf(ind1) = min(temp1, temp2);
                        res.sup(ind1) = max(temp1, temp2);
                    end

                    ind1 = exponent >= 0;
                    temp1 = base.inf(ind1) .^ exponent(ind1);
                    temp2 = base.sup(ind1) .^ exponent(ind1);
                    res.inf(ind1) = min(temp1, temp2);
                    res.sup(ind1) = max(temp1, temp2);

                    ind1 = base.inf < 0 & abs(round(exponent)-exponent) > eps;
                    res.inf(ind1) = NaN;
                    res.sup(ind1) = NaN;
                    
                    ind1 = abs(round(exponent)-exponent) <= eps & base.inf < 0 & ...
                           base.sup > 0 & ~(mod(exponent,2)) & exponent ~= 0;
                    res.inf(ind1) = 0;                    
                end
                
            else
                throw(CORAerror('CORA:dimensionMismatch',base,exponent));
            end            
        end    
    end
    

% a number .^ an inetrval    
elseif isnumeric(base)
    
    res = exponent;
    
    res.inf = min(base .^ exponent.inf, base .^ exponent.sup);
    res.sup = max(base .^ exponent.inf, base .^ exponent.sup);
        
    ind1 = base < 0 ;
    res.inf(ind1) = NaN;
    res.sup(ind1) = NaN;

% an interval .^ an interval   
else
    res = base;
        
        % to find the minimum and the maximum value from all combinations
        possibleValues = {base.inf .^ exponent.inf, base.inf .^ exponent.sup,...
            base.sup .^ exponent.inf, base.sup .^ exponent.sup};
        
        res.inf = possibleValues{1};
        res.sup = possibleValues{1};
        
        for i = 2:4
            res.inf = min(res.inf, possibleValues{i});
            res.sup = max(res.sup, possibleValues{i});
        end
    
    ind1 = base.inf < 0;
    res.inf(ind1) = NaN;
    res.sup(ind1) = NaN;
    
    
end

% return error if NaN occures
if any(any(isnan(res.inf))) || any(any(isnan(res.sup)))
    throw(CORAerror('CORA:outOfDomain','validDomain','base >= 0'));
end

%------------- END OF CODE --------------