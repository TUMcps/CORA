function res = plus(factor1, factor2)
% plus - Overloaded '+' operator
%
% Syntax:
%    res = plus(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - zoo objects
%    order - the cut-off order of the Taylor series. The constant term is
%            the zero order.
%
% Outputs:
%    res - a zoo object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, inteval

% Authors:       Dmitry Grebenyuk
% Written:       05-November-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    if isa(factor1, 'zoo') && isa(factor2, 'zoo')
        
        res = arrayfun(@(a, b) aux_s_plus_zz(a, b), factor1, factor2, 'UniformOutput', false);

    elseif isa(factor1,'zoo') && isa(factor2,'double')

        res = arrayfun(@(a) aux_s_plus_zd(a, factor2), factor1, 'UniformOutput', false);
    
    elseif isa(factor1,'double') && isa(factor2,'zoo')
        
        res = arrayfun(@(b) aux_s_plus_dz(factor1, b), factor2, 'UniformOutput', false);
        
    elseif isa(factor1,'zoo') && isa(factor2,'interval')
        
        res = arrayfun(@(a) aux_s_plus_zd(a, factor2), factor1, 'UniformOutput', false);
        
    elseif isa(factor1,'interval') && isa(factor2,'zoo') 
        
        res = arrayfun(@(b) aux_s_plus_dz(factor1, b), factor2, 'UniformOutput', false);
        
    else
         throw(CORAerror('CORA:wrongValue','first/second',...
             "('zoo','zoo'), ('zoo','double'), ('double','zoo'), ('zoo','interval'), or ('interval','zoo')"));       
        
    end
    A = cat(1, res{:});
    res = reshape(A, size(res));
    
end


% Auxiliary functions -----------------------------------------------------

% Implementation for a scalar

% addition of two zoo objects
function res = aux_s_plus_zz(factor1, factor2)
    [factor1,factor2] = combineZooObjects(factor1,factor2);
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} + factor2.objects{i}; 
    end    
end

% addition of zoo and double
function res = aux_s_plus_zd(factor1, factor2)
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} + factor2; 
    end     
end

% addition of double and zoo
function res = aux_s_plus_dz(factor1, factor2)
    res = factor2;
    for i = 1:length(res.method)
       res.objects{i} = factor2.objects{i} + factor1; 
    end      
end

% ------------------------------ END OF CODE ------------------------------
