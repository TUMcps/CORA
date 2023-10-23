function [CN, state] = batchCombinator(N, K, batch_size, state)
% batchCombinator - Function which generates a sampling without
%    repetitions/replacement of the set 1:N, taken K at a time.
%    Can be called batch-wise to enable large output.
%    Based on the combinator by Matt Fig.
%
% Syntax:
%    [CN, state] = batchCombinator(N, K, batch_size, state)
%
% Inputs:
%    N - sample set size
%    K - sample drawn at a time. Requirement: K<=M
%    batch_size - number of combinations to be returned
%    state - struct with combinator state
%
% Outputs:
%    CN - generated combinations
%    state - new combinator state
%
% Example: 
%    comb_state = struct;
%    while true
%       [new_batch, comb_state] = batchCombinator(30, 20, 1000, comb_state);
%       if comb_state.done == true
%           break;
%       end
%    end
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Michael Eichelbeck
% Written:       27-July-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if ~isfield(state, 'BC') % init
        
        if K>N
            state.done = true;
            CN = [];
            warning("batchCombinator - K>N results in an empty batch")
            return;
        end

        state.M = double(N);  % Single will give us trouble on indexing.
        
        state.WV = 1:K;  % Working vector.
        state.lim = K;   % Sets the limit for working index.
        state.inc = 1;   % Controls which element of WV is being worked on.
        state.stp = 0;   % internal tracker
        state.flg = 0;   % internal tracker
        
        state.BC = floor(prod(state.M-K+1:state.M) / (prod(1:K)));  % total number of combinations
        state.ii = 1; % global index
    
        state.done = false; % have all combinations been returned?

    end 
 
    % copy struct fields to tmp variables for speed
    M = state.M;
    WV = state.WV;
    lim = state.lim;
    inc = state.inc;
    stp = state.stp;
    flg = state.flg;
    BC = state.BC;
    ii = state.ii;
    done = state.done;

    n = batch_size;

    if ii + n > BC % last batch
        n = BC-(ii-1);
        n = floor(n);
    end

    CN = zeros(n,K,class(N));
    
    for nii = 1:n
    
        if ii == BC % The last row. 
            CN(nii,:) = (N-K+1):N;
            done = true;
        elseif ii == 1  % The first row. 
            CN(nii,:) = WV;             
        
        else % all other rows

            if logical((inc+lim)-N) % The logical is nec. for class single(?)
                stp = inc;  % This is where the for loop below stops.
                flg = 0;  % Used for resetting inc.
            else
                stp = 1;
                flg = 1;
            end
            
            for jj = 1:stp
                WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment.
            end
            
            CN(nii,:) = WV;  % Make assignment.
            inc = inc*flg + 1;  % Increment the counter.
            lim = WV(K - inc + 1 );  % lim for next run.     
        
        end % end if
        
        ii = ii + 1;

    end % end for

    % copy local variables to state struct
    state.M = M;
    state.WV = WV;
    state.lim = lim;
    state.inc = inc;
    state.stp = stp;
    state.flg = flg;
    state.BC = BC;
    state.ii = ii;
    state.done = done;
          
end % end of function

% ------------------------------ END OF CODE ------------------------------
