function surf_jet
%SURF_JET Set up extreme color values using COLORMAP JET
%
%Syntax
% surf_jet

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

[X,Y,Z] = sphere(2);
surf(X*0,Y*0,Z*0,-ones(size(Z)),'FaceAlpha',0)
hold on
surf(X*0,Y*0,Z*0,ones(size(Z)),'FaceAlpha',0)
hold on
colormap jet
axis equal
grid off
axis off
if ~haslight(gca)
    camlight right
end
%
%end function