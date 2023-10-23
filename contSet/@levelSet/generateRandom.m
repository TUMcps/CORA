function ls = generateRandom(varargin)
% generateRandom - Generates a random halfspace
%
% Syntax:
%    ls = levelSet.generateRandom()
%    ls = levelSet.generateRandom('Dimension',n)
%
% Inputs:
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       <'NrEquations',nrEqs> - number of equations
%       <'CompOps',compOps> - cell array of operators ('==' or '<=' or '<')
%       <'MaxExponent',maxExp> - maximum exponent in equations
%
% Outputs:
%    ls - random levelSet
%
% Example: 
%    ls1 = levelSet.generateRandom();
%    ls2 = levelSet.generateRandom('Dimension',3);
%    ls = levelSet.generateRandom('Dimension', 2,'NrEquations',2, ...
%           'CompOps',{'<='})
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% name-value pairs -> number of input arguments is always a multiple of 2
if mod(nargin,2) ~= 0
    throw(CORAerror('CORA:evenNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(1:end);
    % check list of name-value pairs
    checkNameValuePairs(NVpairs,{'Dimension','NrEquations','CompOps','MaxExponent'});
    % dimension given?
    [NVpairs,n] = readNameValuePair(NVpairs,'Dimension');
    % Nr of equations given?
    [NVpairs,nrEqs] = readNameValuePair(NVpairs,'NrEquations');
    % possible compare operators given?
    [NVpairs,compOps] = readNameValuePair(NVpairs,'CompOps',{},{'<=','<'});
    % max exponents given?
    [NVpairs,maxExp] = readNameValuePair(NVpairs,'MaxExponent',{},3);
end

% check input arguments
inputArgsCheck({ ...
    {n,'att','numeric','positive'}; ...
    {nrEqs,'att','numeric','positive'}; ...
    {compOps,'att','cell'}; ...
    {maxExp,'att','numeric','positive'}; ...
    });

% default computation for dimension
if isempty(n)
    nmax = 9;
    n = randi(nmax) + 1;
end

% default computation for dimension
if isempty(nrEqs)
    nmax = 4;
    nrEqs = randi(nmax) + 1;
end

% check compOps
if any(cellfun(@(x) ~ismember(x,{'<=','<','=='}),compOps,'UniformOutput',true))
    throw(CORAerror('CORA:wrongValue',"name-value pair 'CompOps'", ...
        "subset of {'==','<=','<'}"));
elseif nrEqs > 1 && any(cellfun(@(x) ismember(x,{'=='}),compOps,'UniformOutput',true))
    throw(CORAerror('CORA:wrongValue',"name-value pair 'CompOps'", ...
         "subset of {'<=','<'} for multiple equations"));
end

% init vars
vars = sym('x_%d',[n 1]);

% generate equations
for i=1:nrEqs
    % A x^e + b <= 0
    A = rand(1,n);
    exps = randi(maxExp,n,1);
    b = randn(1);

    eq(i) = A * vars.^exps + b;
    eqCompOp{i} = compOps{randi(length(compOps))};
end

% extract from cell array for single equation levelSet
if nrEqs == 1
    eqCompOp = eqCompOp{1};
end

% instantiate interval
ls = levelSet(eq,vars,eqCompOp);

% ------------------------------ END OF CODE ------------------------------
