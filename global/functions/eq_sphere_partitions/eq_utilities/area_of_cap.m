function area = area_of_cap(dim, s_cap)
%AREA_OF_CAP Area of spherical cap
%
%Syntax
% area = area_of_cap(dim, s_cap);
%
%Description
% AREA = AREA_OF_CAP(dim, S_CAP) sets AREA to be the area of an S^dim spherical
% cap of spherical radius S_CAP.
%
% The argument dim must be a positive integer.
% The argument S_CAP must be a real number or an array of real numbers.
% The result AREA will be an array of the same size as S_CAP.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, and for dim==3 (when pi/6 <= s_cap <= pi*5/6),
% AREA is calculated in closed form, using the analytic solution of
% the definite integral given in the reference.
% Otherwise, AREA is calculated using the Matlab function BETAINC,
% the incomplete Beta function ratio.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > a=area_of_cap(2,pi/2)
% a =
%     6.2832
% > a=area_of_cap(3,0:pi/4:pi)
% a =
%          0    1.7932    9.8696   17.9460   19.7392
%
%See also
% SRADIUS_OF_CAP

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-24 $ PL
% Use incomplete Beta function BETAINC for dim == 3,
% (when s_cap < pi/6 or s_cap > pi*5/6) and for all dim > 3.
% Use sin(s_cap).^2 in preference to (1-cos(s_cap))/2.
% $Revision 1.01 $ $Date 2005-03-16 $ PL
% Use incomplete Beta function BETAINC for dim > 8.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    area = 2 * s_cap;
case 2
    area = 4*pi * sin(s_cap/2).^2;
case 3
     %
     % Flatten s_cap into a row vector.
     %
     shape = size(s_cap);
     n = prod(shape);
     s_cap = reshape(s_cap,1,n);
     area = zeros(size(s_cap));
     %
     % Near the poles, use the incomplete Beta function ratio.
     %
     pole = (s_cap < pi/6) | (s_cap > pi*5/6);
     area(pole) = area_of_sphere(dim) * betainc(sin(s_cap(pole)/2).^2,dim/2,dim/2);
     %
     % In the tropics, use closed solution to integral.
     %
     trop = s_cap(~pole);
     area(~pole) = (2*trop-sin(2*trop))*pi;

     area = reshape(area,shape);
otherwise
     area = area_of_sphere(dim) * betainc(sin(s_cap/2).^2,dim/2,dim/2);
end
%
% end function
