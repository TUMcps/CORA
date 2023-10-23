function res = times(factor1,factor2)
% times - Overload '.*' operator
%
% Syntax:
%    res = times(factor1,factor2)
%
% Inputs:
%    factor1 - first zoo object
%    factor2 - second zoo object
%
% Outputs:
%    res - resulting zoo object
%
% Other m-files required: interval
% Subfunctions: multiply
% MAT-files required: none
%
% See also: taylm, interval

% Authors:       Dmitry Grebenyuk
% Written:       06-November-2017
% Last update:   ---  
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    if isa(factor1,'zoo') && isa(factor2,'zoo')
            
        res = arrayfun(@(a, b) aux_s_times_zz(a, b), factor1, factor2, 'UniformOutput', false);
    
    elseif isa(factor1,'zoo') && (isa(factor2,'double') || isa(factor2,'interval'))
        
        res = arrayfun(@(a) aux_s_times_zd(a, factor2), factor1,  'UniformOutput', false);
        
    elseif (isa(factor1,'double') || isa(factor1,'interval')) && isa(factor2,'zoo') 

        res = arrayfun(@(b) aux_s_times_dz(factor1, b), factor2, 'UniformOutput', false);
        
    else
         throw(CORAerror('CORA:wrongValue','first/second',...
             "('zoo','zoo'), ('zoo','double'), ('double','zoo'), ('zoo','interval'), ('double','zoo'), or ('interval','zoo')"));
        
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end


% Auxiliary functions -----------------------------------------------------

% Implementation for a scalar
function res = aux_s_times_zz(factor1, factor2)

    [factor1,factor2] = combineZooObjects(factor1,factor2);
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} * factor2.objects{i}; 
    end      
end

function res = aux_s_times_zd(factor1, factor2)

    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} * factor2; 
    end    
end

function res = aux_s_times_dz(factor1, factor2)

    res = factor2;
    for i = 1:length(res.method)
       res.objects{i} = factor2.objects{i} * factor1; 
    end  
end

% ------------------------------ END OF CODE ------------------------------
