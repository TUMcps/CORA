function res = intervalMultiplication(obj,I)
% intervalMultiplication - computes the multiplication of an interval with
% a constraint zonotope. This function is called by the interval/mtimes
% function
%
% Syntax:  
%    res = intervalMultiplication(obj,I)
%
% Inputs:
%    obj - conZonotope object 
%    I - interval object
%
% Outputs:
%    res - conZonotope after multiplication of an interval with a 
%         conZonotope
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    I = interval([0 1; 1 0], [1 2; 2 1]);
%    cZonoRes = I * cZono;
%
%    hold on
%    plot(cZono,[1,2],'r');
%    plot(cZonoRes,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, interval/mtimes, zonotope/intervalMultiplication

% Author:       Niklas Kochdumper
% Written:      03-July-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


    % center and radius of interval matrix
    m=mid(I);
    r=rad(I);

    % absolute value of zontope center and generators
    Zabssum=sum(abs(obj.Z),2);

    % construct resulting conZonotope object
    res = obj;
    res.Z = [m*obj.Z, diag(r*Zabssum)];
    res.A = [obj.A, zeros(size(obj.A,1),size(Zabssum,1))];

    res.ksi = [];
    res.R = [];
 

%------------- END OF CODE --------------