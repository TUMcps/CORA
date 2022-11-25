function res = norm(obj, varargin)
% norm - computes exactly the maximum norm value of all possible matrices
%
% Syntax:  
%    res = norm(obj, varargin)
%
% Inputs:
%    obj - interval matrix
%    varargin - list of optional inputs
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
% See also: @zonotope/norm

% Author:       Matthias Althoff, Victor Gassmann
% Written:      02-November-2017
% Last update:  23-July-2020
% Last revision:---

%------------- BEGIN CODE --------------
type = 2;
if length(varargin)>1
    error('Too many input arguments');
end
if ~isempty(varargin) 
    type = varargin{1};
end
if isnumeric(type) && type==2
    error('Euclidean matrix norm not implemented (exact computation NP hard)');
end

% for type=1,>2; and type='fro' (absolute, see [1]), the matrix attaining
% the maximum respective norm is given by
A_max = abs(center(obj)) + rad(obj);
res = norm(A_max,type);

%------------- END OF CODE --------------