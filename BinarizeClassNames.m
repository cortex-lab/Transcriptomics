function [BinaryNames, ClassRoots] = BinarizeClassNames(Names)
% [BinaryNames, ClassRoots] = BinarizeClassNames(Names)
%
% given a binary class hierarchy specified with periods, e.g.
% {'', '.Cnr1', '.Cnr1', '.Nxph1.Lamp5', '.Nxph1.Sst'}
% turn it into a set of binary strings, e.g.
% {'', '0', '01', '10', '11'}
%
% ClassRoots is a nLevels by 1 cell array, element n is a 2^n cell array of
% strings giving the root for that binary level, or '_' if nothing left
% there

nLevels = max(cellfun(@(x) sum(x=='.'), Names));


nNames = length(Names);

ClassRoots = cell(nLevels,1);
BinaryNames = cell(nNames,1);

Remain = Names;

for Level=1:nLevels+1 % go plus one so they are all given classes on last iter
    ClassRoots{Level} = cellstr(repmat('_',2^Level,1));
    
    % find the next token for all strings
    [ClassToken Remain] = strtok(Remain,'.');
    
    % now go through all roots from the previous level
    for ParentClass=1:2^(Level-1);
        if Level==1
            My = 1:nNames;
            ParentRoot = '';
        else
            ParentRoot = ClassRoots{Level-1}{ParentClass};
            My = strmatch(ParentRoot, Names);
        end
        
        % ParentRoot is the name of the ParentClass
        % My is an index to all those names at this level matching the ParentClass
        
        if isempty(My), continue; end;
        
        % find new tokens at current level
        UniqueTokens = unique(ClassToken(My));
        MyTokens = UniqueTokens(~cellfun(@isempty,UniqueTokens)); %not empty
                
        if length(MyTokens)==2
            % binary class is (ParentClass-1)*2 plus 1 or 2, because MATLAB starts at 1
            ClassRoots{Level}{(ParentClass-1)*2 + 1}=[ParentRoot '.' MyTokens{1}];
            ClassRoots{Level}{(ParentClass-1)*2 + 2}=[ParentRoot '.' MyTokens{2}];
        elseif length(MyTokens)==1
            ClassRoots{Level}{(ParentClass-1)*2 + 1}=[ParentRoot '.' MyTokens{1}];
        elseif length(MyTokens)~=0
            error(sprintf('Neither 0 nor 1 nor 2 tokens under class %s', ClassRoots{Level-1}{ParentClass}));
        end
        
        % now set the class of anyone who doesn't have a token to the
        % parent's name
        NoToken = My(cellfun(@isempty, ClassToken(My)));
        BinaryNames(NoToken) = {dec2bin(ParentClass-1,Level-1)};

        end
    end

end
            