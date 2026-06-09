function [X0, spec, info] = vnnlib2cora(file)
% vnnlib2cora - import specifications from .vnnlib files
%    Dispatches to a version-specific backend based on the
%    (vnnlib-version <X.Y>) header. Files without a header are treated as
%    VNN-LIB 1.x.
%
% Syntax:
%    [X0,spec] = vnnlib2cora(file)
%    [X0,spec,info] = vnnlib2cora(file)
%
% Inputs:
%    file - path to a .vnnlib file storing the specification
%
% Outputs:
%    X0 - cell array of input sets (interval for single-network specs,
%         polytope on the joint input space for multi-network specs)
%    spec - specification object on the (joint) output space
%    info - struct with metadata:
%             .version    - '1.0' or '2.0'
%             .networks   - struct array describing each declared network:
%                             .name, .inputs(j).{name,dtype,shape},
%                             .output.{name,dtype,shape}, .isomorphicTo
%             .isomorphism - struct array of {target, source} pairs
%
% Reference:
%    - https://www.vnnlib.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper, Tobias Ladner, Benedikt Kellner
% Written:       23-November-2021
% Last update:   19-April-2026 (BK, version-dispatched backend)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

version = aux_detectVersion(file);

switch version
    case '1.0'
        [X0, spec] = priv_vnnlib2cora_v1(file);
        if nargout >= 3
            info = aux_v1Info(X0, spec);
        end
    case '2.0'
        [X0, spec, info] = priv_vnnlib2cora_v2(file);
    otherwise
        throw(CORAerror('CORA:notSupported', ...
            sprintf('VNN-LIB version %s is not supported.', version)));
end

end


% Auxiliary functions -----------------------------------------------------

function version = aux_detectVersion(file)
% Scan the file for a (vnnlib-version <X.Y>) header. Strip ;-comments
% first so a header inside a comment is ignored. Default to '1.0'.

text = fileread(file);
text = regexprep(text, ';[^\r\n]*', '');
m = regexp(text, '\(\s*vnnlib-version\s*<\s*([\d\.]+)\s*>\s*\)', ...
    'tokens', 'once');
if isempty(m)
    version = '1.0';
else
    version = m{1};
end

end

function info = aux_v1Info(X0, spec)
% Build a v1-shaped info struct so callers asking for the 3-output form
% on a 1.x file get a consistent contract.

info = struct;
info.version = '1.0';

if isempty(X0)
    inputDims = 0;
else
    inputDims = dim(X0{1});
end
if isempty(spec)
    outputDims = 0;
else
    outputDims = dim(spec(1).set);
end

net = struct;
net.name = '';
net.inputs = struct('name', 'X', 'dtype', 'real', 'shape', inputDims);
net.output = struct('name', 'Y', 'dtype', 'real', 'shape', outputDims);
net.isomorphicTo = '';
info.networks = net;
info.isomorphism = struct('target', {}, 'source', {});

end

% ------------------------------ END OF CODE ------------------------------
