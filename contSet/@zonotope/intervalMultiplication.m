function Z = intervalMultiplication(Z,I)
% intervalMultiplication - computes the multiplication of an interval with
% a zonotope; requires a seperate function since the zonotope class is no 
% longer superior to the interval class
%
% Syntax:  
%    Z = intervalMultiplication(Z,I)
%
% Inputs:
%    Z - zonotope object 
%    I - interval object
%
% Outputs:
%    Z - Zonotpe after multiplication of an interval with a zonotope
%
% Example: 
%    Z=zonotope([1 1 0; 0 0 1]);
%    I=interval([0 1; 1 0], [1 2; 2 1]);
%    plot(Z);
%    hold on
%    Z=intervalMultiplication(Z,I);
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      26-July-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%get center of interval matrix
T=mid(I);
%get radius of interval matrix
S=rad(I);
%auxiliary value
Zabssum=sum(abs(Z.Z),2);
%compute new zonotope
Z.Z=[T*Z.Z,diag(S*Zabssum)]; 
 

%------------- END OF CODE --------------