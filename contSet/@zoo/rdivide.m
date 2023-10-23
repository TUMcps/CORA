function res = rdivide(numerator, denominator)
% rdivide - Overload './' operator
%
% Syntax:
%    res = rdivide(numerator,denominator)
%
% Inputs:
%    numerator - numerator (class zoo, double)
%    denominator - denominator (class zoo, double)
%
% Outputs:
%    res - resulting zoo object
%
% Other m-files required: inverse
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Authors:       Dmitry Grebenyuk
% Written:       06-November-2017
% Last update:   ---  
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    if isa(numerator,'zoo') && isa(denominator,'zoo')
        
        res = arrayfun(@(a, b) aux_s_rdivide_zz(a, b), numerator, denominator, 'UniformOutput', false);
        
    elseif isa(numerator,'zoo') && isa(denominator,'double')
        
        res = arrayfun(@(a) aux_s_rdivide_zd(a, denominator), numerator, 'UniformOutput', false);
        
    elseif isa(denominator,'zoo') && isa(numerator,'double')

        res = arrayfun(@(b) aux_s_rdivide_dz(numerator, b), denominator, 'UniformOutput', false);
        
    end
    
    A = cat(1, res{:});
    res = reshape(A, size(res));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_s_rdivide_zz(numerator, denominator)       
        
    [numerator,denominator] = combineZooObjects(numerator,denominator);
    res = numerator;
    for i = 1:length(res.method)
       res.objects{i} = numerator.objects{i} ./ denominator.objects{i}; 
    end      
end

function res = aux_s_rdivide_zd(numerator, denominator)
        
    res = numerator;
    for i = 1:length(res.method)
       res.objects{i} = numerator.objects{i} ./ denominator; 
    end           
end

function res = aux_s_rdivide_dz(numerator, denominator)

    res = denominator;
    for i = 1:length(res.method)
       res.objects{i} = numerator ./ denominator.objects{i}; 
    end       
end

% ------------------------------ END OF CODE ------------------------------
