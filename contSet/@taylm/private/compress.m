function [res,rest] = compress(obj)
% compress - remove redundant as well as very small coefficients from the 
%            talyor model
%
% Syntax:  
%    [res, rest] = compress(obj)
%
% Inputs:
%    obj - taylm object
%
% Outputs:
%    res - resulting taylm (remainder not updated yet)
%    rest - remainder term resulting from the dicarded coefficients (class
%           interval)
%
% Other m-files required: interval, interval
% Subfunctions: addCoefficients, removeCoefficients
% MAT-files required: none
%
% See also: taylm

% Author:       Niklas Kochdumper, Dmitry Grebenyuk
% Written:      14-June-2017
% Last update:  01-December-2017 (DG)  
% Last revision:---

%------------- BEGIN CODE -------------

    % extract coefficients 
    c = obj.coefficients;
    e = obj.monomials;

    % sort coefficients
    [ e, ind ] = sortrows ( e );
    c = c(ind);
    
    % add coefficients with the same exponent.
    [c,e] = addCoefficients(c,e);
    
    % remove small coefficients and the last coefficients if the order of
    % coefficients is larger then the threshold
    [c,e,cRem,eRem] = removeCoefficients(c,e,obj.max_order,obj.tolerance);
    
    % compute interval overapproximation for the discarded terms
    if isempty(cRem)
        rest = interval(0,0);
    else
        rem = obj;
        rem.coefficients = cRem;
        rem.monomials = eRem;
        rem.remainder = interval(0,0);

        rest = interval(rem);
    end
    
    % remove empty variables
    [e, names_of_var] = removeZeroVar(e, obj.names_of_var);
    
    % assamble resulting taylm
    res = obj;
    res.coefficients = c;
    res.monomials = e;
    res.names_of_var = names_of_var;
    
end


%% Auxiliary functions
    
function [c,e] = addCoefficients(cOld,eOld)

    % initialize output variables
    [len, width] = size(eOld);
    
    c = zeros(len,1);
    e = zeros(len,width);

    % add up coefficients with same exponent
    get = 0;
    put = 0;

    while ( get < len)

        get = get + 1;

        if ( 0 == put )

          put = put + 1;
          c(put) = cOld(get);
          e(put,:) = eOld(put,:);

        else

          if ( all( e(put, :) == eOld(get, :) ) )
            c(put) = c(put) + cOld(get);
          else
            put = put + 1;
            c(put) = cOld(get);
            e(put,:) = eOld(get,:);
          end
        end
    end

    % cut off empty coefficients  
    len = put;
      
    c = c(1:len);
    e = e(1:len,:);
     
end 

function [c,e,cRem,eRem] = removeCoefficients(c,e,N,tolerance)
    
     % intialize remainder coefficients
     [len, width] = size(e);
     
     cRem = zeros(len,1);
     eRem = zeros(len,width);

     % remove small coefficients
     get = 0;
     put = 0;
     putRem = 0;

     while ( get < len )

        get = get + 1;
        
        % N is a max_order
        if ( tolerance < abs ( c(get) ) || e(get, 1) == 0) && e(get, 1) <= N
          put = put + 1;
          c(put) = c(get);
          e(put,:) = e(get,:);
        else
          putRem = putRem + 1;
          cRem(putRem) = c(get);
          eRem(putRem,:) = e(get,:);
        end

     end

     % cut off zeros coefficients
     c = c(1:put);
     e = e(1:put,:);
     
     cRem = cRem(1:putRem);
     eRem = eRem(1:putRem,:);

end

function [e, names] = removeZeroVar(e, names)
    % find zero columns
    ind = any(e,1);
    if any(ind == 0)
        % prune variables
        names = {names{ ind(2:end) }};
        e( :, ~ind ) = [];  %columns
        
        if isempty(names) % to avoid empty list; a dummy variable
            names = {'dummy'};
            e = [0, 0];
        end
    end
end
%------------ END OF CODE ------------
