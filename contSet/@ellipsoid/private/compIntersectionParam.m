function val = compIntersectionParam(W1,q1,W2,q2)
% compIntersectionParam - computes zero root of 'rootfnc' and returns the
%    corresponding argument
%
% Syntax:
%    val = compIntersectionParam(W1,q1,W2,q2)
%
% Inputs:
%    W1,q,W2,q2 - center and shape matrices of E1,E2 (see 'and')
%
% Outputs:
%    val - solution minimizing the volume of the parametrization in 'and'
%          (for details, see [1],[2])
%
% References:
%   [1] Largely based on the Ellipsoidal Toolbox, see:
%       https://de.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
%   [2] For detailed explanations, see:
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and

% Authors:       Victor Gassmann
% Written:       14-October-2019
% Last update:   20-May-2022 (VG, use interval functionality of fzero)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

f = @(p) rootfnc(p,W1,q1,W2,q2);
try
    val = fzero(f,[0,1]);
catch
    % choose the one extremum which has smallest volume
    if det(W1)>det(W2)
        val = 1;
    else
        val = 0;
    end
end


% n = length(q1);
% 
% % since f can be expressed as a polynomial, we build the polynomial by
% % sampling N points (where N is its degree), and then compute its roots,
% % selecting the appropriate (in [0,1]) one
% % f is given by
% %%% y = a*detX^2*trace(X_inv*(W1-W2))-n*detX^2*(2*q'*(W1*q1-W2*q2)+q'*(W2-W1)*q-q1'*W1*q1+q2'*W2*q2)
% % Its degree is given by 2*n-1
% 
% % The interpolation matrix is invertible if and only if we use 4*n
% % distinct sampling points (Vandermonde matrix)
% N = 2*n-1;
% 
% if false
%     X = zeros(N+1,1);
%     for k=0:N
%         X(k+1) = f(exp(-2*pi*1i*k/(N+1)));
%     end
%     c = zeros(N+1,1);
%     for k=0:N
%         c(k+1) = 1/(N+1)*sum(X.*exp(-2*pi*k*1i/(N+1)*(0:N)'));
%     end
% 
% elseif false%true
%     % avoid 0 and 1 due to either W1 or W2 (see e.g. "andHalfspace.m") being
%     % possibly degenerate
%     S = linspace(0.1,0.9,N+1)';
%     
%     % construct interpolation matrix
%     M = S.^(N:-1:0);
%     
%     % find output for all samples (could probably be vectorized)
%     Y = zeros(N+1,1);
%     for k=1:N+1
%         Y(k) = f(S(k));
%     end
%     
%     % find coefficients of polynomial
%     c = M\Y;
% end
% 
% % compute all roots
% r = roots(c);
% 
% % filter out all unwanted roots
% r(~withinTol(abs(imag(r)),0,TOL)) = [];
% r = real(r);
% 
% % remove outside range
% ind = (r<-TOL) | (r>1+TOL);
% if all(ind)
%     % no value in range [0,1] => use closest
%     [~,ii] = min(min([abs(r),abs(r-1)],[],1));
%     val = ii-1;
% else
%     r(ind) = [];
%     % as per ET paper, only one is left
%     % if more than 1 left => choose most "inside" point
%     [~,ii] = max(min([abs(r),abs(r-1)],[],2));
%     val = r(ii);
% end

% ------------------------------ END OF CODE ------------------------------
