function Z = intervalMultiplication(Z,I)
% intervalMultiplication - computes the multiplication of an interval with
%    a zonotope
%
% Syntax:  
%    Z = intervalMultiplication(Z,I)
%
% Inputs:
%    Z - zonotope object 
%    I - interval object
%
% Outputs:
%    Z - zonotope after multiplication of an interval with a zonotope
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    I = interval([0 1; 1 0], [1 2; 2 1]);
%    IZ = intervalMultiplication(Z,I);
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(IZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      26-July-2016 
% Last update:  20-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

%get center of interval matrix
T=center(I);
%get radius of interval matrix
S=rad(I);
%auxiliary value
Zabssum=sum(abs(Z.Z),2);

%compute new zonotope
if ~any(T)
    % no empty generators if interval matrix is symmetric
    Z.Z = [0*center(Z),diag(S*Zabssum)];
else
    Z.Z = [T*Z.Z,diag(S*Zabssum)];
end

%------------- END OF CODE --------------