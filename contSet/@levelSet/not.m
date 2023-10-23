function ls = not(ls)
% not - overloads '~' operator to compute the complement of a level set
%    note: only supported for comparison operators '<' and '<=' and
%    single-equation level set
%
% Syntax:
%    ls = ~ls;
%    ls = not(ls);
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    ls - complement of levelSet object
%
% Example: 
%    % init symbolic variables
%    syms x y z
% 
%    % init level set
%    eq = -x^2 - y^2 + 5;
%    ls = levelSet(eq,[x;y],'<=');
% 
%    % compute complement
%    ls_ = ~ls;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only supported for single-equation level sets
compOp = ls.compOp;
if iscell(compOp)
    if length(compOp) > 1
        throw(CORAerror('CORA:notSupported',['Complement operation is only '...
            'supported for single-equation level sets.']));
    else
        compOp = compOp{1};
    end
end
if ~any(strcmp(compOp,{'<','<='}))
    throw(CORAerror('CORA:notSupported',['Complement operation is only '...
        'supported for comparison operators ''<'' and ''<=''.']));
end

% negate equation
eq_ = -ls.eq;

% change comparisons operators: < to <= and <= to <
if strcmp(compOp,'<=')
    compOp_ = '<';
elseif strcmp(compOp,'<')
    compOp_ = '<=';
end

% init complement
ls = levelSet(eq_,ls.vars,compOp_);

% ------------------------------ END OF CODE ------------------------------
