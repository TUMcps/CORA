function [E] = mtimes(A,E_in)
% mtimes - Overloaded '*' operator for the multiplication of a matrix with
% an ellipsoid
%
% Syntax:  
%    [E] = mtimes(A,E)
%
% Inputs:
%    A - numerical matrix
%    E - Ellipsoid object 
%
% Outputs:
%    E - Ellipsoid 
%
% Example: 
%    E=ellipsoid([1 0; 0 1]);
%    matrix=[1 1; 1 1];
%    plot(E);
%    hold on
%    E=A*E;
%    plot(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Victor Gassmann
% Written:      13-March-2019 
% Last update:  15-October-2019
% Last revision:---

%------------- BEGIN CODE --------------

%Find an ellipsoid object
%Is factor1 an ellipsoid?
if ~isa(E_in,'ellipsoid')
    error('Second factor has to be an ellipsoid object'); 
elseif size(A,2)~=length(E_in.Q)
     if all(size(A)==1)
         A = A*eye(length(E_in.Q));
     else
        error('A must have same number of columns as Q has rows');
     end
end
M = A*E_in.Q*A';
% make sure result is symmetric and eigenvalues are real
E=ellipsoid(1/2*(M+M'),A*E_in.q);
%------------- END OF CODE --------------