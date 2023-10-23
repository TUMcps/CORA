function linSys = generateRandom(varargin)
% generateRandom - Generates a random linear system of the form
%    x' = Ax + Bu
%    y  = Cx
%
% Syntax:
%    linSys = linearSys.generateRandom()
%    linSys = linearSys.generateRandom('Dimension',n)
%    linSys = linearSys.generateRandom('Dimension',n,...
%               'InputDimension',nrInputs,'OutputDimension',nrOutputs)
%    linSys = linearSys.generateRandom('Dimension',n,'RealInterval',realInt)
%    linSys = linearSys.generateRandom('Dimension',n,'RealInterval',realInt,...
%               'ImaginaryInterval',imagInt)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'StateDimension',n> - state dimension
%       <'InputDimension',nrInputs> - input dimension
%       <'OutputDimension',nrOutputs> - output dimension
%       <'RealInterval',realInt> - interval for real part of eigenvalues
%       <'ImaginaryInterval',imagInt> - interval for imaginary part of eigenvalues
%
% Outputs:
%    linSys - random linear system
%
% Example: 
%    linSys1 = linearSys.generateRandom();
%    linSys2 = linearSys.generateRandom('StateDimension',3);
%    linSys3 = linearSys.generateRandom('StateDimension',5,...
%       'RealInterval',interval(-5,-1),'ImaginaryInterval',interval(-1,1));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       07-November-2022
% Last update:   11-November-2022 (MW, integrate into linearSys)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'StateDimension','InputDimension',...
        'OutputDimension','RealInterval','ImaginaryInterval'});
    % state dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'StateDimension');
    % input dimension given?
    [NVpairs,nrInputs] = readNameValuePair(NVpairs,'InputDimension');
    % output dimension given?
    [NVpairs,nrOutputs] = readNameValuePair(NVpairs,'OutputDimension');
    % interval for real part of eigenvalues given?
    [NVpairs,realInt] = readNameValuePair(NVpairs,'RealInterval');
    % interval for imaginary part of eigenvalues given?
    [NVpairs,imagInt] = readNameValuePair(NVpairs,'ImaginaryInterval');
end

% check input arguments or assign default values randomly if not provided
if isempty(n)
    n = randi([4,10]);
end

if isempty(nrInputs)
    nrInputs = randi([1,3]);
end

if isempty(nrOutputs)
    nrOutputs = randi([1,2]);
end

if representsa_(realInt,'emptySet',eps)
    realInt = interval(-1-10*rand(1),-rand(1));
end

if representsa_(imagInt,'emptySet',eps)
    % only one value needed...
    imagMax = 10*rand(1);
    imagInt = interval(-imagMax,imagMax);
end

% check input arguments
inputArgsCheck({{n,'att','numeric','nonnan'};
                {nrInputs,'att','numeric','nonnan'};
                {nrOutputs,'att','numeric','nonnan'};
                {realInt,'att','interval','nonnan'}; ...
                {imagInt,'att','interval','nonnan'}});

% compute value with minimum distance to zero (imaginary parts are always
% complex conjugate)
imagMax = min(abs(infimum(imagInt)),supremum(imagInt));


% check if n is even; if not, at least one eigenvalues HAS TO BE real
nRem = n;
nReal = 0;
if mod(n,2) ~= 0
    nReal = 1;
    nRem = nRem - 1;
    % ...nRem is now even
end

% determine number of complex and real eigenvalues
nConj_max = nRem/2;
nConj = randi([0,nConj_max]);
% resulting remaining eigenvalues are real
nReal = nReal + nRem - 2*nConj;

% we need nReal values
vals_real = unifrnd(infimum(realInt),supremum(realInt),nReal,1);
% additionally, we need nConj real values
val_imag_real = unifrnd(infimum(realInt),supremum(realInt),nConj,1);

% we need nConj imag rand values, other half is complex conjugate
vals_imag = unifrnd(0,imagMax,nConj,1);

% for now, we assume that we have no double eigenvalues (probability of
% this happening is zero anyway)
J_real = diag(vals_real);
J_conj_c = cell(nConj,1);
for i=1:nConj
    J_conj_c{i} = [val_imag_real(i),-vals_imag(i);
                   vals_imag(i),val_imag_real(i)];
end
J = blkdiag(J_real,blkdiag(J_conj_c{:}));

% instantiate random matrix for A​
P = randn(n);
A = P*(J/P);

% input matrix
B = randn(n,nrInputs);

% output matrix
C = randn(nrOutputs,n);

% instantiate linear system​
linSys = linearSys(A,B,[],C);

% ------------------------------ END OF CODE ------------------------------
