function s_cap = sradius_of_cap(dim, area)
%SRADIUS_OF_CAP Spherical radius of spherical cap of given area
%
%Syntax
% s_cap = sradius_of_cap(dim, area);
%
%Description
% S_CAP = SRADIUS_OF_CAP(dim, AREA) sets S_CAP to be the spherical radius of 
% an S^dim spherical cap of area AREA.
%
% The argument dim must be a positive integer.
% The argument AREA must be a real number or an array of real numbers.
% The result S_CAP will be an array of the same size as AREA.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, S_CAP is calculated in closed form.
% Otherwise, S_CAP is approximated using the Matlab function FZERO.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > s_cap=sradius_of_cap(2,area_of_sphere(2)/2)
%
% s_cap =
%     1.5708
%
% > s_cap=sradius_of_cap(3,(0:4)*area_of_sphere(3)/4)
% s_cap =
%          0    1.1549    1.5708    1.9867    3.1416
%
%See also
% FZERO, AREA_OF_CAP

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-25 $
% Use asin rather than acos to avoid subtraction.
% Use symmetry to avoid loss of accuracy in the Southern hemisphere.
% Remove check for when area is close to area of sphere.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    s_cap = area/2;
case 2
    s_cap = 2*asin(sqrt(area/pi)/2);
otherwise
    %
    % Flatten area into a row vector.
    %
    shape = size(area);
    n = prod(shape);
    area = reshape(area,1,n);
    s_cap = zeros(size(area));
    for k = 1:n
        ak = area(k);
        as = area_of_sphere(dim);
        if ak >= as
            s_cap(k) = pi;
        else
            if (2*ak > as)
                ak = as - ak;
                flipped = true;
            else
                flipped = false;
            end
            sk = ...
                fzero(inline(sprintf('area_of_cap(%d,s)-%21.14g',dim,ak),'s'),[0,pi]);
            if flipped
                s_cap(k) = pi - sk;
            else
                s_cap(k) = sk;
            end
        end
    end
    %
    % Reshape output to same array size as original area.
    %
    s_cap = reshape(s_cap,shape);
end
% end function
