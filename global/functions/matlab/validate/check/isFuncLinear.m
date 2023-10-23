function res = isFuncLinear(f,varargin)
% isFuncLinear - checks which equations of a (vector- or matrix-valued)
%    function handle depend linearly on the input arguments
%    caution: this function is only probabilistically correct!
%
% Syntax:
%    res = isFuncLinear(f)
%    res = isFuncLinear(f,inputArgs)
%
% Inputs:
%    f - function handle
%    inputArgs - (optional) vector with size of each input argument to f
%
% Outputs:
%    res - logical array whether (i,j)-th function is a linear function
%
% Example: 
%    f = @(x,u) [x(1) - u(2); x(3)^2*x(2) - u(1)];
%    res = isFuncLinear(f,[3;2])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: withinTol, inputArgsLength

% Authors:       Mark Wetzlinger
% Written:       05-July-2022
% Last update:   06-June-2023 (MW, add quick check for constant functions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if vector length for each input argument of f is given
if nargin == 1
    % read out length of each input argument
    [sizeInputArgs,sizeFunOut] = inputArgsLength(f);
elseif nargin == 2
    sizeInputArgs = varargin{1};
    sizeFunOut = [];
elseif nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% init symbolic input arguments
x = sym('x',[max(sizeInputArgs),length(sizeInputArgs)]);

% quick check: if the function f is constant (which may occur in a
% computation of derivatives), f has one input argument (actually none)
% and the insertion f(x) yields a numeric matrix, no symbolic variables
if isscalar(sizeInputArgs) && sizeInputArgs == 1 && isnumeric(f(x))
    res = true(sizeFunOut); return
end

% check if Hessian matrix is all-zero

% init cell array for evaluation of f
inputArgs = cell(length(sizeInputArgs),1);
% gather all used symbolic variables in one vector for Jacobian evaluation
inputArgsAll = sym('a',[sum(sizeInputArgs),1]);
idxStart = 1;
for i=1:length(inputArgs)
    inputArgs{i} = x(1:sizeInputArgs(i),i);
    inputArgsAll(idxStart:idxStart+sizeInputArgs(i)-1,1) = inputArgs{i};
    idxStart = idxStart + sizeInputArgs(i);
end

% evaluate function to determine number of output arguments and initialize
% resulting logical vector for linear equations
if isempty(sizeFunOut)
    res = false(size(f(inputArgs{:})));
else
    if isscalar(sizeFunOut)
        res = false(sizeFunOut,1);
    else
        res = false(sizeFunOut);
    end
end

% init symbolic all-zero Hessian matrix for comparison
allZeroH = sym(zeros(length(inputArgsAll)));

% loop over each individual function
for j=1:size(res,2)
    % eval (i,j)-th function
    temp = f(inputArgs{:});

    % compute Jacobian
    J = jacobian(temp(:,j),inputArgsAll);
    
    % loop over each equation
    for i=1:size(J,1)
        H = jacobian(J(i,:),inputArgsAll);
        res(i,j) = isequal(H,allZeroH);
    end
end


% % old version based on randomized evaluation
% 
% % linearity of a function depends on two criteria:
% % 1. forall x,y \in domain:                 f(x+y) = f(x) + f(y)
% % 2. for x \in domain and a scalar alpha:   f(alpha*x) = alpha*f(x)
% 
% numInputArgs = length(sizeInputArgs);
% depth = 4;
% 
% % define scalar and vectors to call the function (truncate random numbers
% % to avoid floating-point errors in evaluation as much as possible)
% 
% % random scalar (non-zero)
% alpha = 0;
% while alpha == 0
%     alpha = round(randn,depth);
% end
% 
% % pre-allocate cell arrays for vector calling the function
% x = cell(numInputArgs,1);
% alphax = cell(numInputArgs,1);
% y = cell(numInputArgs,1);
% xy = cell(numInputArgs,1);
% 
% % init counter
% i = 0;
% while i < numInputArgs
%     i = i + 1;
% 
%     % random values for i-th input argument
%     x{i,1} = round(randn(sizeInputArgs(i),1),depth);
%     y{i,1} = round(randn(sizeInputArgs(i),1),depth);
%     % ensure that no values add up to 0: otherwise try other random values
%     if ~all(x{i} + y{i})
%         i = i - 1;
%         continue;
%     end
% 
%     % compute scalar multiplication and addition
%     alphax{i,1} = alpha * x{i};
%     xy{i,1} = x{i} + y{i};
% end
% 
% % evaluate function
% fx = f(x{:});
% fy = f(y{:});
% fxy = f(xy{:});
% falphax = f(alphax{:});
% 
% % check criteria: use relative floating-point accuracy
% res = withinTol(fx+fy,fxy,10*eps) & withinTol(alpha*fx,falphax,10*eps);

% ------------------------------ END OF CODE ------------------------------
