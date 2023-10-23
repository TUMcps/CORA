function displayMatrixVector(x,varname)
% displayMatrixVector - displays a matrix or vector on the command window
%    up to a certain maximum size and abbreviated when all-zero
%
% Syntax:
%    displayMatrixVector(x,varname)
%
% Inputs:
%    x - variable to be displayed (numeric, interval, intervalMatrix)
%    varname - name of variable
%
% Outputs:
%    ---
%
% Example:
%    M = [2 0; -1 2];
%    displayMatrixVector(M,'M');
%
% Other m-files required: DISPLAYDIM_MAX
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/display, linearSysDT/display

% Authors:       Mark Wetzlinger
% Written:       19-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% x ... matrix/vector to be displayed

% get size and text for display
if isnumeric(x)
    x_size = size(x);
    if any(x_size == 1)
        text = "vector";
    else
        text = "matrix";
    end
elseif isa(x,'interval')
    x_size = size(x);
    if any(x_size == 1)
        text = "interval";
    else
        text = "interval matrix";
    end
elseif isa(x,'intervalMatrix')
    x_size = dim(x);
    text = "interval matrix";
elseif isa(x,'matZonotope')
    x_size = dim(x);
    text = 'matrix zonotope';
end

% all-zero check
allzero = false;
if isnumeric(x)
    % standard matrices/vectors
    allzero = ~any(any(x));
elseif isa(x,'interval')
    % interval: all radii have to be zero
    allzero = all(rad(x) == 0);
elseif isa(x,'intervalMatrix')
    % interval matrix: radius has to be zero
    allzero = all(all(rad(x) == 0));
elseif isa(x,'matZonotope')
    % matrix zonotope: all generator matrices have to be all-zero
    allzero = cellfun(@(G) ~any(any(G)),x.generator,'UniformOutput',true);
end

if isempty(x)
    % empty case
    fprintf(newline);
    disp(varname + " = []");
    fprintf(newline);
elseif allzero
    % abbreviation for all-zero matrix/vector
    fprintf(newline);
    disp(varname + " = all-zero " + x_size(1) + "-by-" + x_size(2) + " " + text);
    fprintf(newline);
elseif ~isa(x,'intervalMatrix') && ~isa(x,'matZonotope') ...
        && any(x_size > 1) && diff(x_size) == 0 && all(all(x - eye(x_size) == 0))
    % abbreviation for identity matrices
    fprintf(newline);
    disp(varname + " = " + x_size(1) + "-by-" + x_size(2) + " identity matrix");
    fprintf(newline);
elseif all(x_size <= DISPLAYDIM_MAX)
    % display using correct variable name
    eval(varname + "= x");
else
    % display 
    fprintf(newline);
    disp(varname + " = " + x_size(1) + "-by-" + x_size(2) + " " + text);
    fprintf(newline);
end

% ------------------------------ END OF CODE ------------------------------
