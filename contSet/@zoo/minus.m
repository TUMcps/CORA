function res = minus(factor1, factor2)
% minus - Overloaded '-' operator
%
% Syntax:
%    res = minus(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - zoo objects
%
% Outputs:
%    res - a zoo object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Dmitry Grebenyuk
% Written:       06-November-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
   
    if isa(factor1, 'zoo') && isa(factor2, 'zoo')
        
        res = arrayfun(@(a, b) aux_s_minus_zz(a, b), factor1, factor2, 'UniformOutput', 0);

    elseif isa(factor1,'zoo') && isa(factor2,'double')
        
        res = arrayfun(@(a) aux_s_minus_zd(a, factor2), factor1, 'UniformOutput', 0);
        
    elseif isa(factor1,'double') && isa(factor2,'zoo')
        
        res = arrayfun(@(b) aux_s_minus_dz(factor1, b), factor2, 'UniformOutput', 0);
     
    elseif isa(factor1,'zoo') && isa(factor2,'interval')
        
        res = arrayfun(@(a) aux_s_minus_zd(a, factor2), factor1, 'UniformOutput', 0);
        
    elseif isa(factor1,'interval') && isa(factor2,'zoo') 
        
        res = arrayfun(@(b) aux_s_minus_dz(factor1, b), factor2, 'UniformOutput', 0);
        
    else
         throw(CORAerror('CORA:wrongValue','first/second',...
             "('zoo','zoo'), ('zoo','double'), ('double','zoo'), ('zoo','interval'), or ('interval','zoo')"));
        
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end


% Auxiliary functions -----------------------------------------------------

% Implementation for a scalar
function res = aux_s_minus_zz(factor1, factor2)

    [factor1,factor2] = combineZooObjects(factor1,factor2);
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} - factor2.objects{i}; 
    end 
        
end

function res = aux_s_minus_zd(factor1, factor2)
        
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} - factor2; 
    end  
 
end

function res = aux_s_minus_dz(factor1, factor2)
        
    res = factor2;
    for i = 1:length(res.method)
       res.objects{i} = -factor2.objects{i} + factor1; 
    end 

end

% ------------------------------ END OF CODE ------------------------------
