function [dA,dB] = lin_error2dAB(R,U,hessian,p,varargin)
% lin_error2dAB - computes the uncertainty interval to be added to system
%    matrix set caused by lagrangian remainder of the linearization
%
% Syntax:
%    [dA,dB] = lin_error2dAB(R,U,hessian,p,varargin)
%
% Inputs:
%    R - ???
%    U - ???
%    hessian - ???
%    p - ???
%    ???
%
% Outputs:
%    dA,dB - deviations for A,B caused by lagrange remainder
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
%   -

% Authors:       Victor Gassmann
% Written:       14-May-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

totalInt = interval(R);
inputInt = interval(U);
dim_x = size(totalInt,1);
dim_u = size(inputInt,1);
dxInt = interval(R+(-p.x));
duInt = interval(U+(-p.u));

if length(varargin)==1
    H = hessian(totalInt,inputInt,varargin{1});
else
    H = hessian(totalInt,inputInt);
end

dx = max(abs(dxInt.inf),abs(dxInt.sup));
du = max(abs(duInt.inf),abs(duInt.sup));
dz = [dx;du];
dA = zeros(dim_x,dim_x);
dB = zeros(dim_x,dim_u);
for i = 1:length(H)
    H_ = abs(H{i});
    H_inf = infimum(H_);
    H_sup = supremum(H_);
    H_ = max(H_inf,H_sup);
    M = 1/2*dz'*H_;
    dA(i,:) = M(:,1:dim_x);
    dB(i,:) = M(:,dim_x+1:end);
end

% ------------------------------ END OF CODE ------------------------------
