function res = plus(factor1, factor2)
% plus - Overloaded '+' operator
%
% Syntax:  
%    res = plus(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - zoo objects
%    order  - the cut-off order of the Taylor series. The constat term is
%    the zero order.
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
% See also: taylm, inteval

% Author:       Dmitry Grebenyuk
% Written:      05-November-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    if isa(factor1, 'zoo') && isa(factor2, 'zoo')
        
        res = arrayfun(@(a, b) s_plus_zz(a, b), factor1, factor2, 'UniformOutput', 0);

    elseif isa(factor1,'zoo') && isa(factor2,'double')

        res = arrayfun(@(a) s_plus_zd(a, factor2), factor1, 'UniformOutput', 0);
    
    elseif isa(factor1,'double') && isa(factor2,'zoo')
        
        res = arrayfun(@(b) s_plus_dz(factor1, b), factor2, 'UniformOutput', 0);
        
    elseif isa(factor1,'zoo') && isa(factor2,'interval')
        
        res = arrayfun(@(a) s_plus_zd(a, factor2), factor1, 'UniformOutput', 0);
        
    elseif isa(factor1,'interval') && isa(factor2,'zoo') 
        
        res = arrayfun(@(b) s_plus_dz(factor1, b), factor2, 'UniformOutput', 0);
        
    else
        
        error('Wrong input')
        
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end

%% --------------- Implementation for a scalar --------------

function res = s_plus_zz(factor1, factor2)
        
        [factor1,factor2] = combineZooObjects(factor1,factor2);
        res = factor1;
        for i = 1:length(res.method)
           res.objects{i} = factor1.objects{i} + factor2.objects{i}; 
        end    
end

function res = s_plus_zd(factor1, factor2)

        res = factor1;
        for i = 1:length(res.method)
           res.objects{i} = factor1.objects{i} + factor2; 
        end     
end
    
function res = s_plus_dz(factor1, factor2)
        
        res = factor2;
        for i = 1:length(res.method)
           res.objects{i} = factor2.objects{i} + factor1; 
        end      
end