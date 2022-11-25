function res = minverse(obj)
% minverse - Compute the matrix inverse for Taylor models
%
% This method is implemented from 
% R. Trinchero, P. Manfredi, T. Ding and I. S. Stievano, "Combined 
% Parametric and Worst Case Circuit Analysis via Taylor Models," in IEEE 
% Transactions on Circuits and Systems I: Regular Papers, vol. 63, no. 7, 
% pp. 1067-1078, July 2016.
%
% Syntax:  
%    res = minverse(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm object
%
% Other m-files required: interval, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm

% Author:       Dmitry Grebenyuk
% Written:      18-November-2017
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------

    [m,n] = size(obj);
    r = {obj.monomials};
    c = {obj.coefficients};
    
    if isempty(r)    % Taylor model without polynomial part
        
        res = obj;
        res.remainder = 1/(obj.remainder);
        
    else                % Taylor model with polynomial part
        
        % find constants
        c_f = cellfun(@(r, c) find_const(r, c), r, c, 'UniformOutput', 0);
        c_f = cat(1, c_f{:});
        
        % return the shape of an input object
        c_f = reshape(c_f, m, n);
        
        if isempty(c_f)
            T = obj;
            c_f = 1/0;
            error('division by zero')
        else    
            T = obj - c_f;
            A = inv(c_f);
        end        
        
        for i = 1:n
            unit = zeros(m,n);
            unit(:,i) = 1;
            B = T.*unit;
            
            tr = 1./(1 + trace(A * B));
            A = A - (ones(m, n).* tr) .* (A * B * A);
        end
        res = A; 
    end
end

function c_f = find_const(r, c)
    if all(r(1,1) == 0)                
        c_f = c(1);
    else
        c_f = c( r(:,1) == 0 );
    end
end
%------------ END OF CODE ------------