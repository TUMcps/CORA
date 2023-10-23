function han = plot(fs,varargin)
% plot - plots a projection of an fullspace
%
% Syntax:
%    han = plot(fs)
%    han = plot(fs,dims)
%    han = plot(fs,dims,type)
%
% Inputs:
%    fs - fullspace object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    fs = fullspace(2);
% 
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(fs,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       03-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[fs,dims,NVpairs] = aux_parseInput(fs,varargin{:});

% 2. preprocess
[I,dims] = aux_preprocess(fs,dims);

% 3. plot
han = plot(I,dims,NVpairs{:});

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [O,dims,NVpairs] = aux_parseInput(O,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{O,'att','fullspace'};
                    {dims,'att','numeric',{'nonnan','vector','positive','integer'}}});
    
    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end
    if max(dims) > dim(O)
        throw(CORAerror("CORA:wrongValue",'second','Specified dimensions must not exceed the dimension of the set.'))
    end
    
    % read additional name-value pairs
    if size(dims) == 1
        NVpairs = readPlotOptions(varargin(2:end),'contour');
    else
        NVpairs = readPlotOptions(varargin(2:end),'fill');
    end
end

function [I,dims] = aux_preprocess(fs,dims)

    % project fullspace
    fs = project(fs,dims);
    dims = 1:length(dims);

    % convert to interval
    I = interval(fs);
end

% ------------------------------ END OF CODE ------------------------------
