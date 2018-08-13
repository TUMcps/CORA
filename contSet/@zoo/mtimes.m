function res = mtimes(factor1, factor2)
% mtimes - Overloaded '*' operator for a Taylor model
%
% Syntax:  
%    res = mtimes(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - a taylm objects
%    order  - the cut-off order of the Taylor series. The constat term is
%    the zero order.
%
% Outputs:
%    res - a taylm object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      20-August-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE -------------
    if isa(factor1, 'zoo') && isa(factor2, 'zoo')
        
        [factor1,factor2] = combineZooObjects(factor1,factor2);
        res = factor1;
        for i = 1:length(res.method)
           res.objects{i} = factor1.objects{i} * factor2.objects{i}; 
        end   

    elseif isa(factor1,'zoo') && (isa(factor2,'double') || isa(factor2,'interval'))

        res = factor1;
        for i = 1:length(res.method)
           res.objects{i} = factor1.objects{i} * factor2; 
        end  
    
    elseif (isa(factor1,'double') || isa(factor1,'interval')) && isa(factor2,'zoo')
        
        res = factor2;
        for i = 1:length(res.method)
           res.objects{i} = factor2.objects{i} * factor1; 
        end  
        
    else
        
        error('Wrong input')
        
    end

end
%------------- END OF CODE --------------
