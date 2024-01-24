function res = compareMatrices(M1,M2,varargin)
% compareMatrices - checks if a given matrix has the same columns as
%    another matrix (up to a given tolerance, sign, and potentially in
%    different order) or whether its columns are a subset of the columns of
%    the other matrix (same conditions); we assume no redundancies
%
% Syntax:
%    compareMatrices(M1,M2)
%    compareMatrices(M1,M2,tol)
%    compareMatrices(M1,M2,tol,flag)
%    compareMatrices(M1,M2,tol,flag,ordered)
%    compareMatrices(M1,M2,tol,flag,ordered,signed)
%
% Inputs:
%    M1 - matrix
%    M2 - matrix
%    tol - (optional) tolerance for numerical comparison
%    flag - (optional) type of comparison
%           'equal': M1 has to be exactly M2
%           'subset': M1 only has to be a subset of M2
%           default: 'equal'
%    ordered - (optional) true/false, whether columns have to be in order
%           default: false
%    signed - (optional) true/false, whether columns are equal up to *-1
%
% Outputs:
%    -
%
% Example:
%    M1 = [2 1; 0 2];
%    M2 = [1 2; 2 0];
%    M3 = [2 2; 0 2];
%    compareMatrices(M1,M2) % true
%    compareMatrices(M1,M3) % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: display

% Authors:       Mark Wetzlinger
% Written:       13-November-2022
% Last update:   22-November-2022 (MW, add subset variant)
%                08-May-2023 (TL, ordered)
%                19-January-2024 (MW, add signed)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize result
res = true;

if isempty(M1) && isempty(M2)
    return
end

% parse input arguments
[tol,type,ordered,signed] = setDefaultValues({eps,'equal',false,true},varargin);

% check input arguments
inputArgsCheck({{M1,'att','numeric',{'nonnan','nonempty','finite'}};
                {M2,'att','numeric',{'nonnan','nonempty','finite'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}; ...
                {type,'str',{'equal','subset'}}; ...
                {ordered,'att','logical'}; ...
                {signed,'att','logical'}});

% check if matrices have same number of rows
if size(M1,1) ~= size(M2,1)
    res = false; return
elseif strcmp(type,'equal') && size(M1,2) ~= size(M2,2)
    % number of columns has to be equal
    res = false; return
elseif strcmp(type,'subset') && size(M1,2) > size(M2,2)
    % number of columns cannot be larger
    res = false; return
end

% speed up computation if matrices have to be equal
if ordered && signed && (strcmp(type,'equal') || all(size(M1) == size(M2)))
    res = all(withinTol(M1,M2,tol),'all');
    return
end

% logical indices which columns have been checked (this is faster than
% eliminating the columns from the matrices)
M2logical = false(size(M2,2),1);
% store index to start searching in M2; used for ordered=true
% has to be larger than index of last found column
jmin = 1;

for i=1:size(M1,2)
    % take i-th column of M1 and see if it is part of M2
    found = false;
    for j=jmin:size(M2,2)
        if ~M2logical(j) && (all(withinTol(M1(:,i),M2(:,j),tol)) ...
                || (~signed && all(withinTol(M1(:,i),-M2(:,j),tol))))
            found = true;
            M2logical(j) = true;
            if ordered
                jmin = j+1;
            end
            break
        end
    end
    % exit if no corresponding column found
    if ~found
        res = false; return
    end
end

if strcmp(type,'equal')
    % all columns have to be found
    res = all(M2logical);
elseif strcmp(type,'subset')
    % not all columns have to be found
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
