classdef heap < handle
% heap - min heap
%
% Syntax:
%    obj = nnHelper.heap(heapEntries)
%
% Inputs:
%    heapEntries - cell-array storing keys or struct with key field
%
% Outputs:
%    obj - generated object
%
% Example:
%    heapEntries = {5, 3, 2, 4, 7};
%    heap = nnHelper.heap(heapEntries);
%    disp(heap.pop().key); % = 2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnLayer

% Authors:       Tobias Ladner
% Written:       24-January-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties
    entries = [];
    lastIdx = 0;
end

methods
    function obj = heap(heapEntries)
        % constructor
        n = length(heapEntries);

        for i = 1:n
            entry = heapEntries{i};
            if ~isa(entry, 'struct')
                entryStruct = struct;
                entryStruct.key = entry;
                entryStruct.i = i;
                entry = entryStruct;
            end

            if isempty(obj.entries)
                obj.entries = repmat(entry, 1, n);
            end

            insert(obj, entry)
        end
    end

    function insert(obj, entry)
        % inserts an entry into the heap
        obj.lastIdx = obj.lastIdx + 1;
        obj.entries(obj.lastIdx) = entry;
        siftUp(obj, obj.lastIdx);
    end


    function minEntry = pop(obj)
        % removes the current min entry from the heap
        minEntry = obj.entries(1);
        maxEntry = obj.entries(obj.lastIdx);
        obj.lastIdx = obj.lastIdx - 1;
        obj.entries(1) = maxEntry;
        siftDown(obj, 1);
    end

    function oldMinEntry = replaceMin(obj, newMinEntry)
        % replaces the current min entry with the given entry
        oldMinEntry = obj.entries(1);
        obj.entries(1) = newMinEntry;
        siftDown(obj, 1); % check invariant
    end

    function minEntry = min(obj)
        % returns the current min entry from the heap
        minEntry = obj.entries(1);
    end

    function res = isempty(obj)
        % checks if heap is empty
        res = obj.lastIdx == 0;
    end
end

methods (Access = private)
    function siftDown(obj, pIdx)
        % sifts the entry at index pIdx down until heap invariant is ok
        while true
            pEntry = obj.entries(pIdx);
            % get children
            c1Idx = 2 * pIdx + 0;
            c2Idx = 2 * pIdx + 1;

            if c1Idx > obj.lastIdx
                % no child entry
                return
            end
            if c2Idx > obj.lastIdx
                % only one child
                c2Idx = c1Idx;
            end

            c1Entry = obj.entries(c1Idx);
            c2Entry = obj.entries(c2Idx);

            % check lexicographic order
            [~, idx] = sortrows([pEntry.key; c1Entry.key; c2Entry.key]);
            if idx(1) == 1
                % heap invariant ok
                return
            else
                if idx(1) == 2
                    cIdx = c1Idx;
                    cEntry = c1Entry;
                else
                    cIdx = c2Idx;
                    cEntry = c2Entry;
                end

                obj.entries(pIdx) = cEntry;
                obj.entries(cIdx) = pEntry;

                pIdx = cIdx;
            end
        end
    end

    function siftUp(obj, cIdx)
        % sifts the entry at index cIdx up until heap invariant is ok
        while true
            if cIdx == 1
                % reached min entry
                return;
            end

            cEntry = obj.entries(cIdx);
            pIdx = floor(cIdx/2);
            pEntry = obj.entries(pIdx);

            % check lexicographic order
            [~, idx] = sortrows([pEntry.key; cEntry.key]);
            if idx(1) == 1
                % heap invariant ok
                return
            else
                obj.entries(pIdx) = cEntry;
                obj.entries(cIdx) = pEntry;

                cIdx = pIdx;
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
