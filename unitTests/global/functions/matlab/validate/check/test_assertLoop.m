function res = test_assertLoop
% test_assertLoop - unit test function for assertLoop
%
% Syntax:
%    res = test_assertLoop()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: assertLoop

% Authors:       Tobias Ladner
% Written:       29-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

for i = 1:2
    seed = randi(1000);
    rng(seed)
    msg ='more info';
    
    % test all options
    assertLoop(true,'',[],i);
    assertLoop(true,msg,seed,i);
    assertLoop(true,msg,[],i);
    assertLoop(true,'',seed,i);
    assertLoop(true,i);

    % should fail
    assertThrowsAs(@assertLoop,'',false,'',[],i);
    assertThrowsAs(@assertLoop,'',false,msg,seed,i);
    assertThrowsAs(@assertLoop,'',false,msg,[],i);
    assertThrowsAs(@assertLoop,'',false,'',seed,i);
    assertThrowsAs(@assertLoop,'CORA:wrongValue',false,[],seed,i);
    assertThrowsAs(@assertLoop,'CORA:wrongValue',1,[],seed,i);
    assertThrowsAs(@assertLoop,'CORA:wrongValue',1,i);


    % test second for loop
    for j = 1:2
    
        % test all options
        assertLoop(true,'',[],i,j);
        assertLoop(true,msg,seed,i,j);
        assertLoop(true,msg,[],i,j);
        assertLoop(true,'',seed,i,j);
        assertLoop(true,i,j);
    
        % should fail
        assertThrowsAs(@assertLoop,'CORA:wrongValue',false,[],[],i,j);
        assertThrowsAs(@assertLoop,'',false,msg,seed,i,j);
        assertThrowsAs(@assertLoop,'',false,msg,[],i,j);
        assertThrowsAs(@assertLoop,'CORA:wrongValue',false,[],seed,i,j);
        assertThrowsAs(@assertLoop,'',false,i,j);

    end
end

end

% ------------------------------ END OF CODE ------------------------------
