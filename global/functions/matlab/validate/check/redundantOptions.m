function redundantOptions(varargin)
% redundantOptions - Check if there are redundant setting present in the 
%   options struct
%
% Syntax:
%    redundantOptions(listRed)
%    redundantOptions(options,validOpts)
%
% Inputs:
%    list - cell-array storing the names of the redundant options
%    options - struct storing the algorithm options
%    validOpts - cell-array storing the names of the valid options
%
% Outputs:
%    ---
%
% Example:
%    options.alg = 'dmd';
%    options.zonotopeOrder = 10;
%
%    redundantOptions(options,{'alg'});
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: inputArgsCheck

% Authors:       Niklas Kochdumper
% Written:       06-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    narginchk(1,2);

    % get list of redundant options
    if nargin == 1
        listRed = varargin{1};
    else
        options = varargin{1};
        validOpts = varargin{2};

        allfields = fields(options);
        redIdx = ~ismember(allfields,validOpts);
        listRed = allfields(redIdx);
    end

    % show warning with all redundant options
    redTxt = ['''',strjoin(listRed,''', '''),''''];

    if ~strcmp(redTxt,"''")
        if VALIDATEOPTIONS_ERRORS
            CORAwarning('CORA:redundant',['Redundant options: ' redTxt]);
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
