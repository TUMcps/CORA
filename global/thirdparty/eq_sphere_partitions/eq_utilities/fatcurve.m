function [X,Y,Z] = fatcurve(c,r)
%FATCURVE Create a parameterized cylindrical surface at radius r from curve c
%
%Syntax
% [X,Y,Z] = fatcurve(c,r);
%
%Description
% [X,Y,Z] = FATCURVE(C, R) sets X, Y and Z to be the coordinates of a
% cylindrical surface at radius R from curve C. This function is intended
% for use with the Matlab function SURF to illustrate curves in R^3.
%
%Examples
% > [X,Y,Z] = fatcurve(c,r);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Flesh out description and examples
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(c,1);
n = size(c,2);
m = 8;
h = 0:1/m:1;
phi = h*2*pi;
X = zeros(n,m+1);
Y = X;
Z = Y;
for k = 1:n-1
    u = c(:,k+1)-c(:,k);
    M = null(u');
    if size(M,2) ~= 2
        fprintf('size(M,2) == %d\n',size(M,2));
        M
        c
	return;
    end
    v = M(:,1);
    w = cross(u,v);
    w = w/norm(w);
    if k > 1
        minindex = 0;
        mindist = 2;
        for j = 1:m
            dist = norm(v - circ(:,j));
            if dist < mindist
                mindist = dist;
                minindex = j;
            end
        end
        offs = phi(minindex);
        circ = v*cos(phi-offs) + w*sin(phi-offs);
    else
        circ = v*cos(phi) + w*sin(phi);
    end
    XYZ = c(:,k)*ones(size(phi)) + r*circ;
    X(k,:) = XYZ(1,:);
    Y(k,:) = XYZ(2,:);
    Z(k,:) = XYZ(3,:);
end    
XYZ = c(:,n)*ones(size(phi)) + r*circ;
X(n,:) = XYZ(1,:);
Y(n,:) = XYZ(2,:);
Z(n,:) = XYZ(3,:);
