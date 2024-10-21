function showCORAtipoftheday()
% showCORAtipoftheday - prints the tip of today in the command window
%    Consider adding this function to your startup file to obtain
%    a tip whenever you start Matlab.
%
% Syntax:
%    showCORAtipoftheday
%
% Inputs:
%    -
%
% Outputs:
%    -
%

% Authors:       Tobias Ladner
% Written:       18-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init random number generator
rs = aux_initRandStream();

% get tips
tips = getCORAtipsoftheday(rs);

% select tip
tip = aux_selectTipOfTheDay(rs, tips);

% display tip
aux_showTipMessage(tip);

end


% Auxiliary functions -----------------------------------------------------

function rs = aux_initRandStream()
     % compute from current day

    % Get the current date
    today = datetime('today');
    
    % Convert the date to a number to use as a seed
    seed = convertTo(today,'posixtime');

    % to not interfer with the outside, generate your own random stream
    rs = RandStream('mt19937ar', 'Seed', seed); 
end

function tip = aux_selectTipOfTheDay(rs, tips)
    % select tip
    tip = tips{rs.randi(numel(tips))};
end

function aux_showTipMessage(tip)

    % format tip
    tip = strrep(tip,'\n','\n  ');
    tip = compose(tip);
    tip = tip{1};

    % format tip of the day message
    fprintf('\n')
    fprintf('<strong>CORA - TIP OF THE DAY</strong>\n')
    fprintf('  %s\n', tip)
    fprintf('\n')
end

% ------------------------------ END OF CODE ------------------------------
