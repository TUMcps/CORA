function E = ellipsoid(C,varargin)
% ellipsoid - converts a capsule to an ellipsoid object
%
% Syntax:
%    E = ellipsoid(C)
%    E = ellipsoid(C,method)
%
% Inputs:
%    C - capsule object
%    method - type of conversion: 'exact' (default), 'outer', 'inner'
%
% Outputs:
%    E - ellipsoid
%
% Example: 
%    C = capsule([1; 1], [0.5; -1], 0.5);
%    Eo = ellipsoid(C,'outer');
%    Ei = ellipsoid(C,'inner');
%
%    figure; hold on;
%    plot(C);
%    plot(Eo); plot(Ei);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default method: 'exact'
method = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{C,'att','capsule'},...
                {method,'str',{'exact','outer','inner'}}});

if any(C.g)
    % if the capsule is not a ball (non-zero generator), exact conversion
    % is impossible
    if strcmp(method,'exact')
        throw(CORAerror('CORA:noExactAlg',C));
    end
else
    % reset method since exact conversion is possible
    method = 'exact';
end

% dimension of capsule
n = dim(C);

if strcmp(method,'exact')
    % C.g is an all-zero vector -> capsule is a ball (or empty)
    if representsa_(C,'emptySet',eps)
        E = ellipsoid.empty(n); return
    end
    % square radius to get correct shape matrix
    E = ellipsoid(C.r^2*eye(n),C.c);

elseif strcmp(method,'outer')
    % outer-approximate union of balls at both ends of the capsule
    
    % construct balls as ellipsoids
    E1 = ellipsoid(C.r^2*eye(n), C.c+C.g);
    E2 = ellipsoid(C.r^2*eye(n), C.c-C.g);

    % outer-approximate union (SDP)
    E = or(E1,E2,'outer');

elseif strcmp(method,'inner')
    % map capsule such that generator becomes [1; 0; ...; 0]
    
    % new basis
    [U,~,~] = svd(C.g);
    C_ = C.c + U'*(C + (-C.c));

    % shape matrix
    Q = diag([(vecnorm(C_.g)+C.r)^2;C.r^2*ones(n-1,1)]);

    % instantiate ellipsoid
    E_ = ellipsoid(Q,zeros(n,1));
    % backtransform
    E = U*E_ + C.c;

end


% ------------------------------ END OF CODE ------------------------------
