function C = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with 
%    a capsule
%
% Syntax:
%    C = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or capsule object
%    factor2 - numerical matrix or capsule object
%
% Outputs:
%    C - capsule after multiplication with a matrix
%
% Example: 
%    C = capsule([1; 1], [0; 1], 0.5);
%    M = [0 1; 1 0];
%    C_ = M*C;
%
%    figure; hold on;
%    plot(C);
%    plot(C_);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       04-March-2019
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find a capsule object
[C,matrix] = findClassArg(factor1,factor2,'capsule');

try

    %numeric matrix
    if isnumeric(matrix)
        % new center
        C.c = matrix*C.c;
        % new generator
        C.g = matrix*C.g;
        % new axes of ellipsoid of transformed ball
        newAxes = eig(matrix*matrix');
        C.r = C.r*max(newAxes);
    
    %something else?
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',factor1,factor2));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if representsa_(C,'emptySet',eps)
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
