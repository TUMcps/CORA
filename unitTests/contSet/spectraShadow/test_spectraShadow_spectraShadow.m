function res = test_spectraShadow_spectraShadow
% test_spectraShadow_spectraShadow - unit test function of spectraShadow
%    (constructor)
%
% Syntax:
%    res = test_spectraShadow_spectraShadow
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

% Authors:       Adrian Kulmburg
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2 dimensional box with radius 3 around point [-1;2]:
A0 = eye(3);
A1 = [0 1 0;1 0 0;0 0 0];
A2 = [0 0 1;0 0 0;1 0 0];

% Only A-matrix
S = spectraShadow([A0,A1,A2]);

% Only A and c
S = spectraShadow([A0,A1,A2],[-1;2]);

% A, c, and G
S = spectraShadow([A0,A1,A2],[-1;2],3*eye(2));

% Initialization through ESumRep
S = spectraShadow({[A0,A1,A2], [eye(3)]});

% instantiate H-representation --------------------------------------------

% wrong initializations
A = [1 0 -1];
A_nonSymmetric = A';
ESumRep_ = {A};
c = [0;1];
c_ = [3; 2; 3; 2;];
G_ = eye(3);

% dimension mismatch
assertThrowsAs(@spectraShadow,'CORA:wrongInputInConstructor',A,c_);
assertThrowsAs(@spectraShadow,'CORA:wrongInputInConstructor',A,c,G_);

% incorrect ESumRep structure
assertThrowsAs(@spectraShadow,'CORA:wrongInputInConstructor',ESumRep_);

% empty argument
assertThrowsAs(@spectraShadow,'CORA:wrongValue',[],c);

% too many arguments
assertThrowsAs(@spectraShadow,'CORA:numInputArgsConstructor',[],A,c,G_,c,c);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
