function E = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or
%    scalar with an ellipsoid
%
% Syntax:
%    E = factor1 * factor2
%    E = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - ellipsoid object, numeric matrix/scalar
%    factor2 - ellipsoid object, numeric scalar
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E = ellipsoid([2.7 -0.2;-0.2 2.4]);
%    M = [1 0.5; 0.5 1];
% 
%    figure; hold on;
%    plot(E,[1,2],'b');
%    plot(M*E,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Victor Gassmann
% Written:       13-March-2019 
% Last update:   15-October-2019
%                07-June-2022 (avoid construction of new ellipsoid object)
%                04-July-2022 (VG, input checks, allow E to be class array)
%                05-October-2024 (MW, remove class array)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try

    % matrix/scalar * ellipsoid
    if isnumeric(factor1)
        % -> factor2 must be an ellipsoid object
        E = factor2;
        % compute auxiliary value for new shape matrix
        M = factor1*E.Q*factor1';
        % make sure it is symmetric
        M = 1/2*(M+M');
        
        E.Q = M;
        E.q = factor1*E.q;
        return
    end
    
    % ellipsoid * scalar
    if isnumeric(factor2)
        % -> factor1 must be an ellipsoid object
        E = factor1;
        % compute auxiliary value for new shape matrix
        M = factor2^2*E.Q;
        % make sure it is symmetric
        M = 1/2*(M+M');
        
        E.Q = M;
        E.q = factor2*E.q;
        return
    end

catch ME
    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);
    rethrow(ME);
end

throw(CORAerror('CORA:noops',factor1,factor2));

% ------------------------------ END OF CODE ------------------------------
