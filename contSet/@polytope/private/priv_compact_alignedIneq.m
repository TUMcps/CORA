function [A,b] = priv_compact_alignedIneq(A,b,tol)
% priv_compact_alignedIneq - removes all redundant aligned constraints
%
% Syntax:
%    [A,b] = priv_compact_alignedIneq(A,b,tol)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    tol - tolerance
%
% Outputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% sort the matrix rows to detect aligned normal vectors
[A,ind] = sortrows(A);
b = b(ind);

% remove all aligned halfspaces
counter = 1;
cTemp = 1;
A_ = zeros(size(A));
b_ = zeros(size(b));

while counter < size(A,1)
   
    % check if two normal vectors are identical
    if sum(abs(A(counter,:)-A(counter+1,:))) < tol
        
        a_ = A(counter,:);
        
        % determine minimum b
        bmin = min(b(counter),b(counter+1));
        counter = counter + 2;
        
        while counter <= size(A,1)
            
            if sum(abs(A(counter,:)-a_)) < tol
               bmin = min(bmin,b(counter));
               counter = counter + 1;
            else
               break;
            end
        end
        
        % store unified normal vector
        A_(cTemp,:) = a_;
        b_(cTemp,:) = bmin;
        cTemp = cTemp + 1;
        
    else
        A_(cTemp,:) = A(counter,:);
        b_(cTemp,:) = b(counter,:);
        cTemp = cTemp + 1;
        counter = counter + 1;
    end  
end

% add last element
if counter == size(A,1)
    A_(cTemp,:) = A(end,:);
    b_(cTemp,:) = b(end);
    cTemp = cTemp + 1;
end

% override constraint set
A = A_(1:cTemp-1,:);
b = b_(1:cTemp-1);

% ------------------------------ END OF CODE ------------------------------
