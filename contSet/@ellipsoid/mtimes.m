function E = mtimes(A,E)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
%    an ellipsoid
%
% Syntax:
%    E = mtimes(A,E)
%
% Inputs:
%    A - numerical matrix
%    E - ellipsoid object 
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% if A is scalar, transform to appropriate diagonal matrix
if isnumeric(A) && numel(A)==1
    A = A*eye(dim(E(1)));
end

% check input arguments
inputArgsCheck({{A,'att','numeric',{'ncols',dim(E(1))}};
                {E,'att','ellipsoid'}});

% empty set: result is empty set
if representsa_(E,'emptySet',eps)
    return
end

% make sure all ellipsoids are of same dimension
if ~all(dim(E)==dim(E(1)))
    throw(CORAerror('CORA:wrongValue','second',...
                        'All ellipsoids need same dimension.'));
end

% make sure A is a matrix
if ~ismatrix(A)
    throw(CORAerror('CORA:wrongValue','first','"A" needs to be a matrix.'));
end

for i=1:numel(E)
    % compute auxiliary value for new shape matrix
    M = A*E(i).Q*A';
    % make sure it is symmetric
    M = 1/2*(M+M');
    
    E(i).Q = M;
    E(i).q = A*E(i).q;
end

% ------------------------------ END OF CODE ------------------------------
