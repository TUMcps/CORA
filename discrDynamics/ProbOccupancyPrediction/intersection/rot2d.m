function pts_out = rot2d(theta, pts_in)
% rot2d rotates a set of points of frame in to frame out
% updated: 19-August-2016, Matthias Althoff
% 
assert(size(pts_in,1)==2);
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
pts_out = R*pts_in;
end


