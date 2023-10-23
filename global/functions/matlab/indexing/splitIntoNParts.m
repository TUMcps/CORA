function arr = splitIntoNParts(K,N)
% splitIntoNParts - split the numbers from 1:K into an array of N parts
%
% Syntax:
%    arr = splitIntoNParts(K,N)
%
% Inputs:
%    K - integer number
%    N - integer number
%
% Outputs:
%    arr - Nx2 array (start and end index of each partition)
%
% Example:
%    K = 50;
%    N = 7;
%    arr = splitIntoNParts(K,N);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{K,'att','numeric',{'scalar','integer','positive'}},...
                {N,'att','numeric',{'scalar','integer','positive'}}})

% cannot split K into N>K parts
if N > K
    throw(CORAerror('CORA:wrongValue','second',...
        'has to be smaller or equal to K'));
end

% quick cases
if N == 1
    arr = [1,K]; return
elseif N == K
    arr = [1:K; 1:K]'; return
end

% init array
arr = zeros(N,2);

% main cases
if withinTol(K/N,round(K/N),1e-10)
    % even partition
    shift = round(K/N);
    arr = [1+(0:N-1)*shift; (1:N)*shift]';
    
else
    % uneven partition
    shift = K/N;
    arr = zeros(N,2);
    for i=1:N
        % start of i-th partition
        if i == 1
            arr(i,1) = 1;
        else
            arr(i,1) = arr(i-1,2)+1;
        end
        % end of i-th partition
        if i < N
            arr(i,2) = floor(i*shift);
        else
            arr(i,2) = K;
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
