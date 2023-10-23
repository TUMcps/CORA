function res = norm(intMat,varargin)
% norm - computes exactly the maximum norm value of all possible matrices
%
% Syntax:
%    res = norm(intMat,varargin)
%
% Inputs:
%    intMat - intervalMatrix object
%    type - (optional) p-norm
%
% Outputs:
%    res - resulting maximum norm value
%
% Example: 
%    ---
%
% References:
%    [1]: Raena Farhadsefat, Ji?r´? Rohn and Taher Lotf: Norms of Interval 
%         Matrices (http://uivtx.cs.cas.cz/~rohn/publist/normlaa.pdf)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/norm

% Authors:       Matthias Althoff, Victor Gassmann
% Written:       02-November-2017
% Last update:   23-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(varargin)>1
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% set default values
type = setDefaultValues({2},varargin);

if isnumeric(type) && type==2
    % Euclidean matrix norm not implemented (exact computation NP hard)
    throw(CORAerror('CORA:wrongValue','second','1, >2, ''fro''.'));
end

% for type=1,>2; and type='fro' (absolute, see [1]), the matrix attaining
% the maximum respective norm is given by
A_max = abs(center(intMat)) + rad(intMat);
res = norm(A_max,type);

% ------------------------------ END OF CODE ------------------------------
