function writeCounterexample(fid, x_, y_, vnnlibInfo)
% writeCounterexample - write a counterexample in VNN-LIB 1.x or 2.0 format
%
% Syntax:
%    writeCounterexample(fid, x_, y_, vnnlibInfo)
%
% Inputs:
%    fid - open file handle
%    x_ - input assignment (column vector)
%    y_ - output assignment (column vector)
%    vnnlibInfo - vnnlib2cora metadata (.version, .networks); [] -> v1 format
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: aux_formatShape
% MAT-files required: none
%
% See also: run_instance, vnnlib2cora

% Authors:       Benedikt Kellner
% Written:       05-June-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(fid, ['sat' newline]);

if isempty(vnnlibInfo) || ~isfield(vnnlibInfo, 'version') ...
        || ~strcmp(vnnlibInfo.version, '2.0')
    % v1 legacy: ((X_0 v) ... (Y_0 v) ...)
    fprintf(fid, '(');
    for j = 1:size(x_, 1)
        fprintf(fid, ['(X_%d %f)' newline], j-1, x_(j));
    end
    for j = 1:size(y_, 1)
        fprintf(fid, ['(Y_%d %f)' newline], j-1, y_(j));
    end
    fprintf(fid, ')');
else
    % v2: per-network "name dtype [shape]" header + values, one per line
    xIdx = 1;
    yIdx = 1;
    for netIdx = 1:length(vnnlibInfo.networks)
        net = vnnlibInfo.networks(netIdx);
        for inIdx = 1:length(net.inputs)
            inp = net.inputs(inIdx);
            nElems = prod(inp.shape);
            fprintf(fid, '%s %s %s\n', inp.name, inp.dtype, ...
                aux_formatShape(inp.shape));
            fprintf(fid, '%g\n', x_(xIdx:xIdx+nElems-1));
            xIdx = xIdx + nElems;
        end
        nOut = prod(net.output.shape);
        fprintf(fid, '%s %s %s\n', net.output.name, net.output.dtype, ...
            aux_formatShape(net.output.shape));
        fprintf(fid, '%g\n', y_(yIdx:yIdx+nOut-1));
        yIdx = yIdx + nOut;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function shapeStr = aux_formatShape(shape)
% format [d1 d2 ...] as "[d1,d2,...]"
shapeStr = ['[' strjoin(arrayfun(@num2str, shape, ...
    'UniformOutput', false), ',') ']'];
end

% ------------------------------ END OF CODE ------------------------------
