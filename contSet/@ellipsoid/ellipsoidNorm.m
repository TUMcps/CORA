function res = ellipsoidNorm(E,p)
% ellipsoidNorm - computes the norm of the point p w.r.t. the
%    ellipsoid-norm induced by the ellipsoid E; this is defined similarly
%    to the zonotope-norm defined in [1, Definition 4].
%
% Syntax:
%    res = ellipsoidNorm(E,p)
%
% Inputs:
%    E - ellipsoid
%    p - nx1-array, with n the dimension of E
%
% Outputs:
%    res - ellipsoid-norm of the point p
%
% Example:
%    q = [0;0];
%    Q = diag([3 2]);
%    E = ellipsoid(Q, q);
%    
%    p = [ 0.369 ; 0.685 ];
%    
%    d = ellipsoidNorm(E, p);
%
%    figure; hold on
%    plot(E,[1,2],'b');
%    plot(d*E,[1,2],'g');
%    plot(p(1),p(2),'rx');
%
% References:
%    [1] A. Kulmburg, M. Althoff. "On the co-NP-Completeness of the
%        Zonotope Containment Problem", European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       06-July-2021
% Last update:   26-July-2021 (VG, check for degenerate ellipsoid)
%                04-July-2022 (VG, more precise input check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {p,'att','double',@(p) all(size(p) == [dim(E),1])}});

% The ellipsoid E is given as {x | (x - q)' * Q^(-1) (x - q) <= 1}...
% ... as a consequence, the norm of a point p w.r.t. E would be given by
% p'*Q^(-1)*p, since we have to ignore the center to get a norm.
% So the only thing we need to be careful about is when E is degenerate,
% then we may need to invert Q by hand.
% In theory, we would then need to compute
% res = sqrt(abs(p'*Qinv*p))
% However, this is not quite stable for most cases. That is why we compute
% instead
% q = Q\p;
% res = sqrt(abs(p'*q))
% for the case where Q is invertible; otherwise, we do something similar
% for the degenerate case
if ~isFullDim(E)

     [U,S,Vt] = svd(E.Q);
     s = diag(S);
     s(s<1e-10) = 0;
     Sinv = diag(1./s);
     
     % We consider the ellipsoid E that has been turned by U', such that
     % the axes of E coincide with the canonical ONB; therefore, we also
     % need to rotate p:
     p = U' * p;
     % For numerical reasons, we need to manually set coordinates of p that
     % are very small to zero
     p(abs(p) < 1e-10) = 0;
     
     q = Sinv * p;
     % Now, we need to be careful; Qinv will contain Inf, and p may contain
     % 0, which is the only case that can make the point contained. Thus,
     % we need to set the values of q to zero in this case
     q(isnan(q)) = 0;
else
     q = E.Q\p;
end

res = sqrt(abs(p'*q));
% The square root is just there to make sure that the resulting function is
% a norm (i.e., scaling p by a factor a should yield a|p|, not a^2|p|.


% We need to do one last check: If res = NaN, something weird happened with
% additions and subtractions of Inf. However, this means that some Inf was
% left, so the point has no chance of being contained. Thus, set it
% manually:
if isnan(res)
    res = Inf;
end

% ------------------------------ END OF CODE ------------------------------
