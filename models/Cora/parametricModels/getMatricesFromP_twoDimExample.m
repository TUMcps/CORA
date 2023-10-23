function [A, B, C, D] = getMatricesFromP_twoDimExample(p,sys)
% getMatricesFromP_twoDimExample - returns matrices of a linear system
% given a parameter vectors. This is required for gray-box conformance
% syntehsis of linear systems.
%
% Syntax:  
%    [A, B, C, D] = getMatricesFromP_twoDimExample(p)
%
% Inputs:
%    p - parameter vector
%    sys - linear sysstem
%
% Outputs:
%    A - system matrix
%    B - input matrix
%    C - output matrix
%    D - matrix
%

% Author:       Matthias Althoff
% Written:      30-August-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain sizes of matrices
dim_x = size(sys.A,1);
dim_u = size(sys.B,2);
dim_y = size(sys.C,1);

% system matrix
A = reshape(p(1:dim_x^2), dim_x, dim_x); 
% remove already used parameters
p = p(dim_x^2+1:end);

% input matrix
B = reshape(p(1:dim_u*dim_x), dim_x, dim_u);
% remove already used parameters
p = p(dim_u*dim_x+1:end);

% output matrix
C = reshape(p(1:dim_x*dim_y), dim_y, dim_x);
% remove already used parameters
p = p(dim_x*dim_y+1:end);

% feedthrough matrix
D = reshape(p(1:dim_u*dim_y), dim_y, dim_u);

end