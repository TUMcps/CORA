function res = mtimes(factor1, factor2)
% mtimes - Overloaded '*' operator for a Taylor model
%
% Syntax:
%    res = mtimes(factor1, factor2)
%
% Inputs:
%    factor1 and factor2 - a taylm objects
%    order - the cut-off order of the Taylor series. The constant term is
%            the zero order.
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

% Authors:       Dmitry Grebenyuk
% Written:       20-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isa(factor1, 'zoo') && isa(factor2, 'zoo')
    % zoo and zoo
    
    [factor1,factor2] = priv_combineZooObjects(factor1,factor2);
    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} * factor2.objects{i}; 
    end   

elseif isa(factor1,'zoo') && (isa(factor2,'double') || isa(factor2,'interval'))
    % zoo and double/interval

    res = factor1;
    for i = 1:length(res.method)
       res.objects{i} = factor1.objects{i} * factor2; 
    end  

elseif (isa(factor1,'double') || isa(factor1,'interval')) && isa(factor2,'zoo')
    % double/interval and zoo
    
    res = factor2;
    for i = 1:length(res.method)
       res.objects{i} = factor2.objects{i} * factor1; 
    end  
    
else % throw error
     throw(CORAerror('CORA:wrongValue','first/second',"'double', 'interval', or 'zoo'"));
end

% ------------------------------ END OF CODE ------------------------------
