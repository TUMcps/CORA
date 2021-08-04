function T = vecalign(x,y)
% computes T such that x || T*y (x parallel to T*y)
[U1,~,V1] = svd(x);
[U2,~,V2] = svd(y);
T = U1*V1*V2'*U2';
end