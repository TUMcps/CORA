function table_authors()
% table_authors - generates a table as a categorical matrix. The rows are the
% name of the authors. The columns are the names of the files. The content
% of each cell is the number of contributions of each author in each file.
% The symbol '-' indicates no contribution. 
%
% Syntax:  
%    table_authors()
%
% Inputs:
%    -
%
% Outputs:
%    -

% Author:       Tobias Ladner
% Written:      20-October-2023
% Last update:  ---
% Last revision:---

% ------------------------------ BEGIN CODE -------------------------------

% extracting the main files in the directory ------------------------------

mainDirs = dir(CORAROOT);
mainDirs = mainDirs([mainDirs.isdir]); %extracting the files that are directories

% remove '.', '..', '.git', and 'manual' directories
unwantedDirs = {'.','..','.git','manual','codegen'};
keepIdx = true(length(mainDirs),1);
for i=1:length(mainDirs)
    keepIdx(i) = ~any(strcmp(mainDirs(i).name,unwantedDirs));
end
mainDirs = mainDirs(keepIdx);

% init --------------------------------------------------------------------

dict = dictionary();

% search authors in all files ----------------------------------------------

disp("Searching for authors:")

% iterate over directories
for i=1:length(mainDirs)
    fprintf(" /%s\n", mainDirs(i).name)
    
    % find *.m files within directory
    files = findfiles([mainDirs(i).folder filesep mainDirs(i).name]);
    files = excludefiles(files, ['global' filesep 'thirdparty']);
    
    % iterate over files
    for j=1:length(files)
         % read content of file
        filetext = fileread([files(j).folder filesep files(j).name]);
    
        % split lines
        lines = splitlines(filetext);

        % find '% Author' line
        for l=1:length(lines)
            line = lines{l};

            if startsWith(line,'% Author')
                % parse authors
                authors = aux_parseAuthors(line);

                % iterate over authors
                for a=1:length(authors)
                    % update dictionary
                    
                    author = authors{a};
                    if ~dict.isConfigured || ~isKey(dict,author)
                        % init author dict
                        dict(author) = aux_initAuthorDict(mainDirs);
                    end

                    % increase count
                    author_dict = dict(author);
                    author_dict(mainDirs(i).name) = author_dict(mainDirs(i).name) + 1;
                    dict(author) = author_dict;
                end

                break
            end
        end
    end
end
disp(' ')

% postprocessing ----------------------------------------------------------

% filter unknown authors
filter_authors = ["???","---"];
authors = keys(dict);
disp("Removing unknown authors")
for a = 1:length(authors)
    if ismember(authors{a},filter_authors)
        fprintf('- ''%s'' (%d)\n',authors{a},sum(values(dict(authors{a}))));
        dict(authors{a}) = [];
    end
end
disp(' ')

% sort by total contribution ----------------------------------------------

authors = keys(dict);
totalConts = zeros(size(authors));

for a=1:length(authors)
    totalConts(a) = sum(values(dict(authors{a})));
end

[totalConts,idx] = sort(totalConts,'descend');
authors = authors(idx);

% print latex table -------------------------------------------------------

threshold = 15;

% init table
fprintf("Latex table for manual:\n\n")
table = CORAtable('latex', ... % switch to 'html' for website
    [{'Name'},arrayfun(@(dir) dir.name,mainDirs,'UniformOutput',false)'], ...
    [{'s'},repmat({'i'},1,numel(mainDirs))], ...
    'ColumnWidths', [max(arrayfun(@strlength,authors)),repmat(5,1,numel(mainDirs))] ...
);
table.printHeader();

% always have Matthias first
aux_printtableline(table,dict, 'Matthias Althoff');

% print remaining authors above threshold
for a=1:length(authors)
    if totalConts(a) >= threshold
        author = authors{a};
        if ~strcmp(author,'Matthias Althoff')
            aux_printtableline(table,dict, author);
        end
    else
        % only list names below threshold
        break;
    end
end
table.printFooter();

% print all others by name
fprintf('\n\nWe also want to thank ')
fprintf('%s, ', authors{a:end-1})
fprintf('and %s for individual contributions.\n\n', authors{end})

disp('.. please update the manual and the website.')

end


% Auxiliary functions -----------------------------------------------------

function aux_printtableline(table,dict,author)
    table.printContentRow([{author},num2cell(values(dict(author)))'])
end

function author_dict = aux_initAuthorDict(mainDirs)

    % init author dict and pre-populate with 0s
    author_dict = dictionary();
    for i=1:length(mainDirs)
        author_dict(mainDirs(i).name) = 0;
    end

end

function authors = aux_parseAuthors(line)
    authors = strtrim(line(12:end));
    authors = strsplit(authors,', ');
end
 
% ------------------------------ END OF CODE ------------------------------
