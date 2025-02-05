function Y = priv_outputSet_canonicalForm(linsys,X,V,v,k)
% priv_outputSet_canonicalForm - computes the output set for linear systems 
%    in canonical form y = Cx + v; we additionally use the step k as an
%    input argument to easily index the vector array v which may be a
%    single vector
%
% Syntax:
%    Y = priv_outputSet_canonicalForm(linsys,X,V,v,k)
%
% Inputs:
%    linsys - linearSys object
%    X - time-point reachable set at time t_k, or
%        time-interval reachable set over time interval [t_k,t_k+1]
%    V - measurement uncertainty
%    v - measurement uncertainty vector
%    k - step
%
% Outputs:
%    Y - output set
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% index vector: if it is an array of vector, we use the step k, otherwise
% it is constant and we can just use it as it is
if size(v,2) > 1
    v = v(:,k);
end

% evaluate y = Cx + v
Y = linsys.C * X + V + v;

% ------------------------------ END OF CODE ------------------------------
