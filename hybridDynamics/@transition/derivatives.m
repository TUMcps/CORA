function trans = derivatives(trans,varargin)
% derivatives - compute derivatives for nonlinear reset functions
%
% Syntax:
%    trans = derivatives(trans)
%    trans = derivatives(trans,fpath)
%    trans = derivatives(trans,fpath,fname)
%
% Inputs:
%    trans - transition object or class array
%    fpath - path to generated file
%    fname - file name
%
% Outputs:
%    trans - updated transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearReset/derivatives

% Authors:       Mark Wetzlinger
% Written:       15-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments, set defaults, check values
narginchk(1,3);
[fpath,fname] = setDefaultValues({...
    [CORAROOT filesep 'models' filesep 'auxiliary' filesep 'transition'],...
    'nonlinear_reset_function'},varargin);
inputArgsCheck({{trans,'att','transition'};...
                {fpath,'att','char'};...
                {fname,'att','char'}});

% loop over transitions
for i=1:numel(trans)
    % update .reset (only if nonlinearReset object)
    if isa(trans(i).reset,'nonlinearReset')
        trans(i).reset = derivatives(trans(i).reset,...
            fpath,sprintf('transition_%i_%s',i,fname));
    end
end

% ------------------------------ END OF CODE ------------------------------
