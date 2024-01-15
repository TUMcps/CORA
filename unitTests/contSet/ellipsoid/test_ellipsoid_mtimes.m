function res = test_ellipsoid_mtimes
% test_ellipsoid_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_ellipsoid_mtimes
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%NOTICE: Before executing this test, make sure test_ellipsoid_supportFunc and
%test_ellipsoid_contains are sucessful as this test makes use of both
%functions.

res = false;
% create ellipsoid
E1 = ellipsoid.generateRandom('IsDegenerate',true);
%If actual dimension is ==1
if E1.dim>1
    %Generate random points within E1
    N = 1000;
    samples = randPoint(E1,N);
    if ~contains(E1,samples)
        disp('test_ellipsoid_mtimes failed');
        return;
    end

    A = randn(ceil(abs(randn)),length(E1.Q));
    %First compute the ellipsoid using framework
    EA = A*E1;
    %Compute the linear transformation using the sample points
    samples_A = A*samples;

    if ~contains(EA,samples_A)
        disp('test_ellipsoid_mtimes failed');
        return;
    end
end
%%%Simple matrix multiply test
Q = [1,0.1,0.8;0.1,3,2;0.8,2,5];
q = [10;-23;15];
E1 = ellipsoid(Q,q);
A = [1,2,3;4,5,6];
EA = A*E1;
true_result=ellipsoid([87.2,193.7;193.7,433.4],[9;15]);

% empty set
% if ~representsa_(mtimes(A,ellipsoid.empty(3)),'emptySet',eps)
%     res = false; return
% end

res = isequal(EA,true_result);

% ------------------------------ END OF CODE ------------------------------
