% GeneSet: a MATLAB object for analysis of Single-Cell RNASeq Data
% (C) Kenneth D. Harris, 2014-2017
% Any questions or comments, please email kenneth.harris@ucl.ac.uk.
%
% This code may rely on some additional library functions. If you get
% "undefined function or variable" errors, please email me and I will sent
% them.
%
% Released under the GPL

classdef GeneSet
    properties
        nGenes; % number of genes
        nCells; % number of cells
        GeneExp; % nGenes x nCells array giving expression levels
        GeneName; % Ngenes x 1 cell array with name of each gene
        GeneInfo; % struct array giving attributes of each gene
        CellInfo; % struct array giving attributes of each cell
        Class; % string giving class assignment of every cell
        CellName; % a name for every cell, to help you remember them.
        Excluded; % a GeneSet containing any 
        tSNE; % a 2 by nCells array giving the output of tSNE dim reduction
        
        GeneCorrels; % nGenes x nGenes array to speed up MDS. Calculated only when needed
    end
    methods

        function g = GeneSet(FileName, nHeaderLines) 
            % g = GeneSet(FileName, nHeaderLines) 
            % load data from CEF or loom file 
            % if passed a cell array of strings, will load all those files
            % if FileName is a .tsv file, will load that treating
            % nHeaderLines as headers
            
            if nargin<1
                return;
            end
            
            % load all in cell array 
            if iscell(FileName)
                g = GeneSet(FileName{1});
                for i=2:length(FileName)
                    
                    g = g.horzcat(GeneSet(FileName{i}));
                end
                % give new names to ensure no duplicates
                g.CellName=RandomNames(g.nCells);
                return
            end

            % find file name extension
            if ~char(FileName)
                error('Was expecting either a filename or cell array of filenames');
            end
            
            LastPeriod = find(FileName=='.', 1, 'last');
            Extension = FileName(LastPeriod:end);

            % there are two file formats: loom and cef
            % can also load a .tsv file if you insist
            if strcmpi(Extension, '.loom')
                fprintf('Loading loom file %s\n', FileName);
                g = GeneSet;

                g.GeneName = deblank(h5read(FileName, '/row_attrs/Gene'));
                % for some reason, some files encase the gene name in
                % "b'GName'"
%                 if g.GeneName{1}(1)=='b'
%                     for i=1:length(g.GeneName)
%                         g.GeneName{i} = g.GeneName{i}(3:end-1);
%                     end
%                 end
                
                g.GeneExp = h5read(FileName, '/matrix')'; % transpose it for nGenes by nCells!
                [g.nGenes, g.nCells] = size(g.GeneExp);

                ColInfo = h5info(FileName, '/col_attrs');
                RowInfo = h5info(FileName, '/row_attrs');

                g.GeneInfo = struct;
                for i=1:length(RowInfo.Datasets);
                    FieldName = RowInfo.Datasets(i).Name;
                    AttrName = FieldName;
                    AttrName(FieldName==' ') = '_'; % get rid of spaces
                    if AttrName(1)=='_'
                        AttrName = AttrName(2:end);
                    end
                    g.GeneInfo = setfield(g.GeneInfo, AttrName, h5read(FileName, ['/row_attrs/' FieldName]));
                end

                g.CellInfo = struct;
                for i=1:length(ColInfo.Datasets);
                    FieldName = ColInfo.Datasets(i).Name;
                    AttrName = FieldName;
                    AttrName(FieldName==' ') = '_'; % get rid of spaces
                    if AttrName(1)=='_'
                        AttrName = AttrName(2:end);
                    end
                    g.CellInfo = setfield(g.CellInfo, AttrName, h5read(FileName, ['/col_attrs/' FieldName]));
                    g.CellInfo.Loom_Row = (1:g.nCells)';
                    g.CellInfo.Loom_File = repmat({FileName},[g.nCells,1]); 
                end

                g.CellName = RandomNames(g.nCells);
            elseif strcmpi(Extension, '.cef')
                fprintf('Loading cef file %s\n', FileName);
                %first line in cef file is info about size of each part
                fid = fopen(FileName);
                FirstLine = textscan(fid,'%s %d %d %d %d %d %d %*[^\r\n]', 1, 'Delimiter', '\t');
                if FirstLine{1}{1}~='CEF'
                    error('CEF file should start with string "CEF"');
                end
                [nHead, nRowAttr, nColAttr, nRows, nCols, Flags] = FirstLine{2:7};
                g.nGenes = nRows;
                g.nCells = nCols;

                % read in header
                Header = textscan(fid, '%s %s %*[^\r\n]', nHead, 'Delimiter', '\t');

                % read in Column Attributes
                g.CellInfo = struct;
                for i=1:nColAttr
                    % read attribute name, skipping multiple blank cells first
    %                 ColAttrName = textscan(fid, '%s', 1, 'Delimiter', '\t', 'MultipleDelimsAsOne', 1);
    %                 ColAttrVals = textscan(fid, '%s', nCols, 'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
    %                 g.CellInfo = setfield(g.CellInfo,ColAttrName{1}{1}, ColAttrVals{1});
                    textscan(fid, '%*[\r\n]'); Line = fgetl(fid); % skip blanks and read line - had to do it this way to avoid some weird file problems
                    [ColAttrName, pos] = textscan(Line, '%s', 1, 'Delimiter', '\t', 'MultipleDelimsAsOne', 1);
                    if pos<length(Line)
                        ColAttrVals = textscan(Line(pos+1:end), '%s', nCols, 'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
                        g.CellInfo = setfield(g.CellInfo,ColAttrName{1}{1}, ColAttrVals{1});
                    end

                end

                % read in Row Attributes
                RowAttrName = textscan(fid, '%s', nRowAttr, 'Delimiter', '\t');
                RowAttrVals = textscan(fid, [repmat('%s ',[1 nRowAttr]) '%*[^\r\n]'], nRows);

                g.GeneName = RowAttrVals{1};
                g.GeneInfo = struct;
                for i=2:nRowAttr
                    g.GeneInfo = setfield(g.GeneInfo, RowAttrName{1}{i}, RowAttrVals{i});
                end

                fclose(fid);
                % now load in actual data
                g.GeneExp = dlmread(FileName, '\t', 3+nColAttr, 1+nRowAttr);
                % strip off blank rows and columns
                g.GeneExp = g.GeneExp(1:g.nGenes, 1:g.nCells);

                % give the cells some memorable names
                g.CellName=RandomNames(g.nCells);
            elseif strcmpi(Extension, '.tsv')
                fprintf('Loading tsv file %s with %d header lines\n', FileName, nHeaderLines);
                
                % read in actual data
                g.GeneExp = dlmread(FileName, '\t', nHeaderLines,1);
                g.nGenes = size(g.GeneExp,1);
                
                
                % read in header lines
                fid = fopen(FileName, 'r');
                for i=1:nHeaderLines
                    Line = fgetl(fid);
                    SplitLine = strsplit(Line, '\t');
                    g.CellInfo = setfield(g.CellInfo, SplitLine{1}, SplitLine(2:end)');
                    if i==1
                        g.nCells = length(SplitLine)-1;
                    elseif g.nCells ~= length(SplitLine)-1
                        error('header lines not all same length')
                    end
                end

                % read gene names
                tmp = textscan(fid, ['%s %*[^\r\n]'], g.nGenes);
                g.GeneName= tmp{1};
                fclose(fid);

                % give them random names
                g.CellName=RandomNames(g.nCells);

                
            else
                error('Unknown file extension');
            end

        end
        
        function IDs = NamesToIDs(g, Names);
            % GeneIDs = NamesToIDs(g, Names);
            % returns numeric IDs of genes whose names are in cell array Names
            % if Names is a string, will just return one matching this
            % if Names is numeric, it will just return it as-is
            % any non-matches are dropped with a warning
            
            if isstr(Names)
                IDs = strmatch(Names, g.GeneName, 'exact'); 
                if isempty(IDs)
                    warning(sprintf('Could not find gene of name %s', Names));
                end
            elseif iscell(Names)
%                 [~, ia, IDs] = intersect(Names, g.GeneName, 'stable');
%                 if length(ia)<length(Names)
%                     str = ['Could not find gene of name ' sprintf('%s ', Names{setdiff(1:length(Names), ia)})];
%                     warning(str);
%                 end
                  [Good, IDs] = ismember(Names, g.GeneName);
                  if any(~Good)
                        str = ['Could not find gene of name ' sprintf('%s ', Names{~Good})];
                        warning(str);
                  end
                
            elseif isnumeric(Names)
                IDs = Names;
            else
                error('Names must be cell or string or numeric');
            end
            
            IDs = IDs(IDs>0);
        end
        
        function IDs = NameStartsWith(g, Root)
            % IDs = NameStartsWith(Root)
            % 
            % finds gene IDs for genes starting with Root.
            % Root can be a string, or a cell array of strings
            
            if isstr(Root)
                IDs = strmatch(Root, g.GeneName);
            elseif iscell(Root)
                IDs = [];
                for i=1:length(Root)
                    IDs = union(IDs, strmatch(Root{i}, g.GeneName));
                end
            else
                error('Root should be string or cell');
            end
        end
            
        function Names = IDsToNames(g, IDs);
            % GeneNames = IDsToNames(g, IDs);
            % returns names of genes with specified numeric IDs 
            % if you pass an array it returns a cell array of strings
            % if you pass a number it returns a string
            
            if isstr(IDs) | iscell(IDs)
                Names = IDs;
            else
                Names = g.GeneName(IDs);
            end
            if iscell(Names) & length(Names)==1
                Names = Names{1};
            end
        end
        
        function IDs = MeanExceeds(g, Threshold)
            % MeanExceeds(g, Threshold)
            % returns IDs of all genes whose mean expression exceeds Threshold
            
            MeanExp = mean(g.GeneExp,2);
            IDs = find(MeanExp>Threshold);
            
        end
        
        function x= Exp(g, Gene, Cells);
            % x = Exp(Gene, Cells);
            % return a vector of expression levels for Gene
            % 
            % Default for Cells is everyone; you can also pass names,
            % classes, etc. (See IdentifyCells);
            % Gene can be numeric, string, or cell array of strings
            % any strings that don't match gene names, will try entry in
            % CellInfo
            
            if nargin<3
                Cells = 1:g.nCells;
            else
                Cells = g.IdentifyCells(Cells);
            end
            
            if isstr(Gene)
                Gene = {Gene};
            end
            
            if isnumeric(Gene)
                x = g.GeneExp(Gene,Cells);
            else
                IsaGene = ismember(Gene, g.GeneName);
                GeneID = g.NamesToIDs(Gene(IsaGene));
                if isstruct(g.CellInfo)
                    IsaCellInfo = ismember(Gene, fieldnames(g.CellInfo));
                else
                    IsaCellInfo = zeros(size(IsaGene));
                end
                
                x = nan(length(Gene), length(Cells));
                x(IsaGene,:) = g.GeneExp(GeneID,Cells);
                if any(IsaCellInfo) % because 0-by-1 matrix runs for loop
                    for i = find(IsaCellInfo)
                        x(i,:) = getfield(g.CellInfo, Gene{i});
                    end
                end
            end
        end
        
        function h = GoodExpressers(g, ThisManyCells, ExpressThisMuch)
            % h = GoodExpressers(g, ThisManyCells, ExpressThisMuch)
            % returns a subset structure with all genes for which at least
            % ThisManyCells cells Express at least ExpressThisMuch
            
            Sorted = sort(g.GeneExp,2, 'descend');
            IDs = find(Sorted(:,ThisManyCells)>=ExpressThisMuch);
            h = g.GeneSubset(IDs);
           
        end
        
        function h = GeneSubset(g, IDs)
            % h = GeneSubset(g, GeneIDs)
            % make a new array corresponding subset of genes
            % if names, they will be converted
            
            IDs = g.NamesToIDs(IDs);
            
            h = g;
            h.GeneExp = g.GeneExp(IDs,:);
            h.GeneName = g.GeneName(IDs);
            h.nGenes = length(IDs);
        end
        
        function h = ExcludeGenes(g, ExcludeThese)
            % h = ExcludeGenes(g, ExcludeThese)
            % 
            % removes the specified genes, and puts them in the Excluded
            % field
            %
            % ExcludeThese can be:
            %   a string, in which case it has to match exactly
            %   a string ending in *, in which case it matches the start
            %   a numeric array, in which case it uses these IDs
            %   a cell array, in which case it takes the union of all
            %   members
            
            if ~iscell(ExcludeThese)
                ExcludeThese = {ExcludeThese};
            end
            
            ExcludeList = [];
            for i=1:length(ExcludeThese)
                e = ExcludeThese{i};
                if isnumeric(e)
                    eIDs = e;
                elseif isstr(e) & e(end)=='*'
                    eIDs = g.NameStartsWith(e(1:end-1));
                elseif isstr(e)
                    eIDs = g.NamesToIDs(e);
                else
                    error('ExcludeGenes: unexpected type');
                end
                
                ExcludeList = union(ExcludeList, eIDs);
            end
            
            h = g.GeneSubset(setdiff(1:g.nGenes, ExcludeList));
            h.Excluded = g.GeneSubset(ExcludeList);
        end
        
        function h = DisexcludeGenes(g)
            % h = DisexcludeGenes(g)
            % puts the excluded genes back in the main structure
            
            h1 = g;
            h1.Excluded = [];
            h = vertcat(h1, g.Excluded);
        end
        
        function h = ReClass(g, NewClass);
            % h = ReClass(g, NewClass);
            % changes the Class attribute to NewClass
            %
            % if a string, uses that structure of CellInfo
            % if numeric, convert to strings (no space padding)
            % if a containers.Map, use that to convert 
            % otherwise uses as-is
            
            h = g;
            if isempty(NewClass)
                h.Class = repmat({''}, g.nCells,1);
            elseif isstr(NewClass)
                h = g.ReClass(getfield(g.CellInfo, NewClass));
            elseif isnumeric(NewClass) | islogical(NewClass)
                h.Class = cellstr(num2str(NewClass(:), '%-d'));
            elseif isa(NewClass, 'containers.Map');
                [uClass, ~, ClassNo] = unique(g.Class, 'stable');
                for i=1:max(ClassNo)
                    h.Class(ClassNo==i) = {NewClass(uClass{i})};
                end
            else
                h.Class = NewClass;
            end
        end
        
        function h = MedianInClass(g)
            % h = MedianInClass(g)
            % replaces expression of every gene with the median in that class
            
            h = g;
            [ClassNames, ~, ClassNum] = unique(g.Class);
            for k=1:length(g.Class);
                MyCells = find(ClassNum==k);
                h.GeneExp(:,MyCells) = repmat(median(g.GeneExp(:,MyCells),2), [1 length(MyCells)]);
            end
        end
            
        
        function [h, Cells] = CellSubset(g, varargin)
            % [h, CellIDs] = CellSubset(g, Cells, Exclude)
            % make a new array corresponding subset of cells
            % see IdentifyCells for more info on how the Cells argument is
            % parsed
            %
            % optional output argument Cells tells you which ones were
            % picked
            
            Cells = g.IdentifyCells(varargin{:});

            
            h = GeneSet;
            h.GeneExp = g.GeneExp(:,Cells);
            h.GeneName = g.GeneName;
            h.nGenes = g.nGenes;
            h.nCells = length(Cells);
            if ~isempty(g.CellName), h.CellName = g.CellName(Cells); end
            if ~isempty(g.tSNE), h.tSNE = g.tSNE(:,Cells); end

            
            if ~isempty(g.CellInfo)
                fn = fieldnames(g.CellInfo);
                for i=1:length(fn)
                    Data = getfield(g.CellInfo, fn{i});
                    h.CellInfo = setfield(h.CellInfo, fn{i}, Data(Cells));
                end
            end
            
            if ~isempty(g.Class)
                h.Class = g.Class(Cells);
            end
            
            if ~isempty(g.Excluded)
                h.Excluded = g.Excluded.CellSubset(Cells);
            end

        end
        
        function h = SortByClass(g, ClassOrder)
            % h = SortByClass(ClassOrder)
            % 
            % reorders cells by their class
            % ClassOrder should be a cell array of strings saying what
            % order the classes should be in. Default is alphabetical
            
            if nargin<2 
                ClassOrder = unique(g.Class);
            end
            [~, ClassNumber] = ismember(g.Class, ClassOrder);
            [~, CellOrder] = sort(ClassNumber);
            
            h = g.CellSubset(CellOrder);
        end
        

        
        function h = Permute(g, Cells)
            % h = Permuate(g, Cells)
            %
            % randomizes expression of all genes, only between the cells in Cells
            % (default is all of them)
            
            % preprocessing to get cell list
            if nargin<2
                Cells = 1:g.nCells;
            else
                Cells = g.IdentifyCells(Cells);
            end
                        
            % now randomize
            h = g;
            nCols = length(Cells);
            for i=1:h.nGenes
                h.GeneExp(i,Cells) = h.GeneExp(i,Cells(randperm(nCols)));
            end
        end
        
        function h = Randomize(g)
            % h = Randomize(g)
            %
            % Randomizes gene expression using Patefield's algorithm, so
            % that all cells and genes have same total expression
            
            h = g;
            % find find and exclude zero genes
            NonZeroGenes = find(any(h.GeneExp>0,2));
            
            % make expression submatrix for these
            Expression = h.GeneExp(NonZeroGenes,:);
            % now apply Patefield's algorithm
            [~,~,RandomizedExpression,Error] = Patefield(size(Expression,1), size(Expression,2), sum(Expression, 2), sum(Expression, 1), 0, 0);            
            if (Error); error('Patefield algorithm error'); end
            h.GeneExp(NonZeroGenes,:) = RandomizedExpression;
        end
            

        function h = Power(g, p);
            % Power(p);
            % transform expression by raising each gene to power p (usually<1)
            h = g;
            h.GeneExp = g.GeneExp.^p;
        end
        
        function h = Log(g, c);
            % Log(c);
            % transform expression by log(x + c)
            
            h = g;
            h.GeneExp = log(g.GeneExp+c);
        end
        
        function h = Clip(g, Val, Genes);
            % Clip(Val, Genes)
            % 
            % maxes out gene expression at Val for genes in Genes (default
            % all)
            
            h = g;
            if nargin<3
                Genes = 1:g.nGenes;
            else
                Genes = g.NamesToIDs(Genes);
            end
            h.GeneExp(Genes,:) = min(h.GeneExp(Genes,:),Val);
        end
        
        function h = Saturate(g, Val, Genes);
            % Saturate(Val, Genes)
            % 
            % applies a Naka-Rushton function to make expression saturate
            % at Val (default 50). Runs only on genes Genes.
            %
            % if you pass Val=inf, it won't do anything.
            
            if nargin<2
                Val = 50; 
            end
            
            if ~isfinite(Val)
                h = g;
                return;
            end
            
            if nargin<3
                Genes = 1:g.nGenes;
            else
                Genes = g.NamesToIDs(Genes);
            end
            
            h = g;

            h.GeneExp(Genes,:) = Val*h.GeneExp(Genes,:)./(Val + h.GeneExp(Genes,:));
        end
        
        function h = ScaleGene(g, p, q);
            % transform expression by dividing the expression of each gene
            % by its norm raised to power p
            % optional input q uses the Lq norm (default 1);
            
            if nargin<3; q=1; end;
                         
            h = g;
            Norm = mean(g.GeneExp.^q,2).^(1/q);
            h.GeneExp = bsxfun(@rdivide, g.GeneExp, Norm.^p + eps);
        end
        
        function h = ScaleCell(g, p, q);
            % transform expression by dividing the expression of each gene
            % by its mean value raised to power p, then rescaling so total
            % expression of entire set is constant.
            % 
            % optional input q uses the Lq norm (default 1);
            
            if nargin<3; q=1; end;

            Norm = mean(g.GeneExp.^q,1).^(1/q);
            Not0Norm = (Norm>0);
            h = g.CellSubset(Not0Norm);

            h.GeneExp = bsxfun(@rdivide, g.GeneExp(:,Not0Norm), Norm(Not0Norm).^p);
            
            ScaleFac = (sum(g.GeneExp(:).^q)/sum(h.GeneExp(:).^q)).^(1/q);
            h.GeneExp = h.GeneExp*ScaleFac;
        end
        
        function h = vertcat(varargin)
            % concatenate the genes in multiple GeneSets into one
            
            h = varargin{1};
            for i=2:length(varargin)
                h.GeneExp = [h.GeneExp; varargin{i}.GeneExp];
                h.GeneName = [h.GeneName; varargin{i}.GeneName];
                h.GeneInfo = [h.GeneInfo; varargin{i}.GeneInfo];
                h.nGenes = h.nGenes + varargin{i}.nGenes;
            end
        end
        
        function h = horzcat(g1, g2, varargin)
            % concatenate the cells in multiple GeneSets into one
            % note that unlike Merge, this DOES require the genes are in the
            % same order (which it checks) But it is a lot faster.
            % 
            % if any CellInfo fields don't match, they are discarded
            
            if nargin>2 % recursion if more than 2 inputs
                h = horzcat(horzcat(g1, g2), varargin{:});
                return;
            end
            
            h = GeneSet;
            
            h.tSNE = []; % this is meaningless when you merge datasets


            if ~isequal(g1.GeneName, g2.GeneName)
                error(sprintf('Gene Names mismatch'));
            end
            h.nGenes = g2.nGenes;
            h.GeneName = g2.GeneName;
            h.GeneInfo= g2.GeneInfo;

            h.GeneExp = [g1.GeneExp, g2.GeneExp];
            h.nCells = g1.nCells + g2.nCells;
            h.Class = vertcat(g1.Class(:), g2.Class(:));
            h.CellName = vertcat(g1.CellName(:), g2.CellName(:));
            
            if ~isempty(g1.CellInfo) & ~isempty(g2.CellInfo)
                % merge cell info (this is the hard part)
                CellInfoFields = intersect(fieldnames(g1.CellInfo), fieldnames(g2.CellInfo));
                NewCellInfo = struct;
                for i=1:length(CellInfoFields(:));
                     f = CellInfoFields{i};
    %                 if isfield(g2.CellInfo, f)
                    c1 = getfield(g1.CellInfo,f);
                    if length(c1)~=g1.nCells
                        error(sprintf('CellInfo.%s has %d entries instead of %d', ...
                            f, length(c1), g1.nCells));
                    end

                    c2 = getfield(g2.CellInfo,f);
                    if length(c2)~=g2.nCells
                        error(sprintf('Arg1 CellInfo.%s has %d entries instead of %d', ...
                            f, length(c2), g2.nCells));
                    end

                    if isnumeric(c1) & ~isnumeric(c2) % this is really annoying
                        c1 = cellstr(num2str(c1));
                    elseif isnumeric(c2) & ~isnumeric(c1) % this is really annoying
                        c2 = cellstr(num2str(c2));
                    end

                    NewCellInfo = setfield(NewCellInfo, f, [c1; c2]);

    %                 else
    %                     fprintf('Deleting CellInfo field %s\n', f);
    % %                     h.CellInfo = rmfield(h.CellInfo, f);
    %                 end
                end

                h.CellInfo = NewCellInfo;
                MissingFields = setdiff(union(fieldnames(h.CellInfo), fieldnames(g2.CellInfo)), CellInfoFields);
                if ~isempty(MissingFields)
                    wString = ['Dropping CellInfo fields ' ...
                        sprintf('%s ', MissingFields{:}) 'which are not in all inputs'];
                    warning(wString);
                end
            end
            
        end
        
        function l = Merge(g1, varargin);
            % h = g1.merge(g2, ...);
            % merge together the cells in two or more GeneSets
            % 
            % note that the genes are NOT assumed to be in the same order
            % use this for example to merge GEO microarray data into GeneSets
            % 
            % for every gene in g2 that matches one in g1 by name, an entry 
            % is added. All other genes set to NaN
            %
            % optional last argument says which genes to keep:
            % 'any' (default): any gene in any input; if not present, set to NaN
            % 'first': keep only genes from g1
            % 'all': keep only genes in all genesets
            
            if any(strcmp(varargin(end), {'any', 'all', 'first'}))
                Option = varargin(end);
                varargin = varargin(1:end-1);
            else
                Option = 'any';
            end
            
            % get list of all genes to keep
            GeneList = g1.GeneName(:);
            for i=1:length(varargin)
                if strcmp(Option, 'any')
                    GeneList = union(GeneList(:), varargin{i}.GeneName);
                elseif strcmp(Option, 'all')
                    GeneList = intersect(GeneList(:), varargin{i}.GeneName);
                end
            end
            
            % now reprocess each one
            gNew{1} = g1;
            gNew(2:length(varargin)+1) = varargin;
            for i=1:length(gNew)
                [MyGenes, MyIndices] = ismember(GeneList, gNew{i}.GeneName);
                NewExp = nan(length(GeneList), gNew{i}.nCells);
                NewExp(MyGenes,:) = gNew{i}.GeneExp(MyIndices(MyGenes),:);
                gNew{i}.GeneExp = NewExp;
                gNew{i}.GeneName = GeneList;
                gNew{i}.nGenes = length(GeneList);
            end
           
            % now use horzcat
            l = horzcat(gNew{:});
        end
                
            
            
            
%             if nargin>=3
%                 l = g1.Merge(g2).Merge(varargin{:});
%                 return;
%             end
%             
%             h = g1;
%             h.nCells = g1.nCells + g2.nCells;
%             h.GeneExp = zeros(h.nGenes, h.nCells);
%             h.GeneExp(:,1:g1.nCells) = g1.GeneExp;
%             h.CellName = vertcat(h.CellName(:), g2.CellName(:));
%             % should copy stuff for CellInfo from horzcat
% 
%             NewRange = (g1.nCells+1):h.nCells;
%             ArrayGenes = zeros(g1.nGenes,1);
%             for i=1:g1.nGenes
%                 Matches = find(strcmp(g1.GeneName{i}, g2.GeneName));
%                 if length(Matches)>=1
%                     h.GeneExp(i,NewRange) = sum(g2.GeneExp(Matches,:),1);
%                     ArrayGenes(i) = 1;
%                 else
%                     h.GeneExp(i,NewRange) = NaN;
%                 end
%             end
%             
%             l = h.GeneSubset(find(ArrayGenes));
            
%         end
        
        
        function h = AddUnitGene(g)
            % add a new gene containing all ones to the bottom of the list
            % (for use in regression etc)
            
            h = g;
            h.nGenes = g.nGenes+1;
            h.GeneExp(h.nGenes,:) = 1;
            h.GeneName{h.nGenes} = 'Unit';
        end
        
        function h = AddGenes(g, Expn, Names)
            % h = AddGenes(g, Expn, Names)
            %
            % add pseudogenes to a Geneset. 
            % Expn: nAddGenes x nCells matrix of expression levels
            % Names: nAddGenes length cell array of gene names (strings)
            
            
            h = g;
            h.GeneExp = [g.GeneExp; Expn];
            h.GeneName = vertcat(g.GeneName, Names(:));
        end
        

           
        function Scatter(g, xGene, yGene, Group, Jitter, varargin);
            % Scatter(xGene, yGene, Group, Jitter, varargin);
            %
            % plot a scatter plot. xGene is a numeric index or name of the
            % gene to plot on the x-axis, yGene same for y axis.
            % or they can just be an array of length nCells.
            %
            % Group determines colors (see gscatter)
            %
            % Jitter defaults to 0.3
            % varargin is passed to gscatter
    
           
            if nargin<5
                Jitter = 0.3;
            end
            
            if isnumeric(xGene) & length(xGene)==g.nCells
                xData = xGene(:)';
                xl = 'Data';
            else
                %xData = g.GeneExp(g.NamesToIDs(xGene),:);
                xData = g.Exp(xGene);
                xl = g.IDsToNames(xGene);
            end
            
            if isnumeric(yGene) & length(yGene)==g.nCells
                yData = yGene(:)';
                yl = 'Data';
            else
                %yData = g.GeneExp(g.NamesToIDs(yGene),:);
                yData = g.Exp(yGene);
                yl = g.IDsToNames(yGene);
            end

            nPoints = length(xData);
            if nargin<4 | isempty(Group)
                if ~isempty(g.Class)
                    Group = g.Class;
                else
                     Group = ones(nPoints, 1);
                end
            end
            
            % symbols: alternate between four of them to help when colors
            % are similar
            nGroups = length(unique(Group));
            ColorMap = HsvNotYellow(ceil(nGroups*1.2)); % 1.2 is to avoid confusion of beginning and end
            Sym = repmat('o+*hxsd^<v>p', [1 ceil(nGroups/12)]);
            
           
            hold off
            
            if isempty(Group) | length(Group)==1
                 gscatter(xData+rand(1,g.nCells)*Jitter, yData+rand(1,g.nCells)*Jitter, varargin);
            else
                gscatter(xData+rand(1,g.nCells)*Jitter, yData+rand(1,g.nCells)*Jitter, Group, ColorMap, Sym', varargin);%, 'interpreter', 'none');
            end
            xlabel(xl, 'Interpreter', 'none');
            ylabel(yl, 'Interpreter', 'none');
            
            
            % now fit line
            if g.nCells>0
                [b bint r rint stats] = regress(yData(:), [xData(:), ones(nPoints,1)]);
                [rho p] = corr(yData(:), xData(:), 'type', 'Spearman');

                title(sprintf('Psn: y = %.2fx + %.0f.  R^2 = %.2f.  p = %.3f. Spmn R=%.2f, p=%.3f.\n', b, stats([1 3]), rho, p));

                hold on;

                plot(xlim, xlim*b(1) + b(2),'k');
            end

        end
        
        function [h, ax, bigax] = PlotMatrix(g, IDs, Group, Jitter, ColorMap, Sym, varargin);
            % [h, ax, bigax] = PlotMatrix(IDs, Group, Jitter, ColorMap, Sym, varargin);
            %
            % makes a scatter plot matrix of genes with specified IDs
            % IDs can be gene numbers or names. If it is a 2-element cell array, will do a
            % scatter of one vs. the other rather than a square matrix
            %
            % Group says what color to plot. If it has length the number of
            % cells, it means one per cell. Otherwise, it should be a cell
            % array specifying a list of class branch names to highlight in color,
            % on top of all others in gray. Defaults to g.Class
            % 
            % Jitter defaults to .3
            % Colors defaults to hsv for the number of classes, 
            % Sym defaults to alternating. If you pass a GeneSet to Colors, it will 
            % choose colors and symbols according to that - so do if you
            % call with a subset
            %
            % varargin is passed to gplotmatrix - so Siz, Doleg, Dispopt, ...
            
            if nargin<3 | isempty(Group)
                if ~isempty(g.Class)
                    Group = g.Class;
                else
                    Group = ones(g.nCells, 1);
                end
            end
            
            if length(Group)<g.nCells
                % assign groups according to input classes
                GpNames = Group;
                % assign all cells by default to group 0
                Group = cell(g.nCells,1);
                Group(1:g.nCells) = {'..'};
                % assign all others
                WhichWhere = zeros(g.nCells, 1);
                for i=1:length(GpNames)
                    MyCells = strncmp(GpNames{i},g.Class, length(GpNames{i}));
                    Group(MyCells)=GpNames(i);
                    WhichWhere(MyCells) = i;
                end 
                % now put them in the order of the classes 
                [~,order] = sort(WhichWhere);
                g = g.CellSubset(order);
                Group = Group(order);
            end
            
            if nargin<4 | isempty(Jitter)
                Jitter = 0;
            end
            
            nGroups = length(unique(Group));
            
            if nargin<5 | isempty(ColorMap);
%                 ColorMap = hsv(nGroups);
%                 TooYellow = find(ColorMap(:,1)>.7 & ColorMap(:,2)>.7);
%                 ColorMap(TooYellow,:) = ColorMap(TooYellow,:)*diag([.5 .5 .7]);
                ColorMap = HsvNotYellow(ceil(nGroups*1.2));
            end

            if nargin<6 | isempty(Sym)
                Sym = repmat('o+svp', [1 ceil(nGroups/6)]);
            end
            
            if strcmp(class(ColorMap), 'GeneSet')
                nG0 = length(unique(ColorMap.Class));
                ColorMap0 = HsvNotYellow(ceil(nG0*1.2));
                Sym0 = repmat('o+svp', [1 ceil(nG0/6)]);
                 
                [~, MyClasses] = ismember(unique(Group, 'stable'), unique(ColorMap.Class,'stable'));
                ColorMap = ColorMap0(MyClasses,:);
                Sym = Sym0(MyClasses);
            end
            
            

            
            if length(IDs)== 2 & iscell(IDs) 
%                 IDs1 = g.NamesToIDs(IDs{1});
%                 IDs2 = g.NamesToIDs(IDs{2});
                IDs1 = IDs{1};
                if isstr(IDs1); IDs1 = {IDs1}; end
                IDs2 = IDs{2};
                if isstr(IDs2); IDs2 = {IDs2}; end
                nIDs1 = length(IDs1);
                nIDs2 = length(IDs2);
                [h, ax, bigax] = gplotmatrix(g.Exp(IDs1)' + rand(g.nCells, nIDs1)*Jitter, ...
                    g.Exp(IDs2)' + rand(g.nCells, nIDs2)*Jitter, Group, ColorMap, Sym, varargin{:});
                set(ax, 'box', 'off');
                for i=1:nIDs1
                    xlabel(ax(nIDs2,i), g.IDsToNames(IDs1(i)), 'Interpreter', 'none'); 
                end
                for i=1:nIDs2
                    ylabel(ax(i,1), g.IDsToNames(IDs2(i)), 'Interpreter', 'none'); 
                end
            elseif size(IDs,2)==1 | size(IDs,1)==1
                % single vector of IDs
%                 IDs = g.NamesToIDs(IDs);
%                 nIDs = length(IDs);
% 
%                 [h, ax, bigax] = gplotmatrix(g.GeneExp(IDs,:)' + rand(g.nCells, nIDs)*Jitter, ...
%                     [], Group, ColorMap, Sym, varargin{:});
                nIDs = length(IDs);

                [h, ax, bigax] = gplotmatrix(g.Exp(IDs)' + rand(g.nCells, nIDs)*Jitter, ...
                    [], Group, ColorMap, Sym, varargin{:});
                set(ax, 'box', 'off');
                for i=1:nIDs
                    xlabel(ax(nIDs,i), g.IDsToNames(IDs(i)), 'Interpreter', 'none'); 
                    ylabel(ax(i,1), g.IDsToNames(IDs(i)), 'Interpreter', 'none'); 
                end
            else
                error('PlotMatrix needs a vector of IDs or a 2-element cell array');
            end

        end
        
        function Ggobi(g, Genes, Group);
            % GGbobi(IDs)
            % launch a ggobi process to visualize specified genes
            
            Genes = g.NamesToIDs(Genes);
            GeneNames = g.IDsToNames(Genes);
            nGenes = length(Genes);
            
            if nargin<3 | isempty(Group)
                if ~isempty(g.Class)
                    Group = g.Class;
                else
                    Group = cellstr(num2str(ones(g.nCells,1)));
                end
            end
            
            if isnumeric(Group) | islogical(Group)
                Group = cellstr(num2str(Group(:)));
            end
            
            if isempty(g.CellName)
                CellName = cellstr(num2str([1:g.nCells]'));
            else
                CellName = g.CellName;
            end
            
            docNode = com.mathworks.xml.XMLUtils.createDocument('ggobidata');
            ggobidata = docNode.getDocumentElement;
            
            %activeColorScheme = docNode.createElement('activeColorScheme');
            %activeColorScheme.setAttribute('name', 'Set3 12');
            %ggobidata.appendChild(activeColorScheme);
            nColors = 9;
            
            
            data = docNode.createElement('data');
            data.setAttribute('name', 'Gene Expression')
            

            variables = docNode.createElement('variables');
            data.appendChild(variables);
            variables.setAttribute('count', num2str(nGenes));

            for v=1:nGenes
                var = docNode.createElement('realvariable');
%                 name = docNode.createElement('name').appendChild(docNode.createTextNode(GeneNames{v}));
%                 var.appendChild(name);
                var.setAttribute('name', GeneNames{v});
                variables.appendChild(var);
            end

            % get unique number of each Cell's Group, which will be its color
            [c ia Color] = unique(Group);
            MaxColor = max(Color);
            
            records = docNode.createElement('records');
            data.appendChild(records);
            records.setAttribute('count', sprintf('%d', g.nCells));
            for r=1:g.nCells
                rec = docNode.createElement('record');
                records.appendChild(rec);
                rec.setAttribute('label', sprintf('%s: %s', CellName{r}, Group{r}));
                rec.setAttribute('color', sprintf('%d', mod(Color(r), nColors)));
%                 if MaxColor>nColors
                    rec.setAttribute('glyphType', sprintf('%d', 5-floor(Color(r)./ nColors)));
                    rec.setAttribute('glyphSize', '1');
%                 end
                rec.appendChild(docNode.createTextNode(sprintf('%f ', g.GeneExp(Genes,r))));
            end

            ggobidata.appendChild(data);
%             xmlwrite(docNode)
            xmlwrite('gobidata.xml', docNode);
             system('start ggobi gobidata.xml');
            
%             fp = fopen('gobidata.csv', 'w');
%             for v=1:length(Names);
%                 fprintf(fp, '"%s"', Names{v});
%                 if v<length(Names)
%                     fprintf(fp, ',');
%                 else
%                     fprintf(fp, '\n');
%                 end
%             end
%             fclose(fp);
%             
%             dlmwrite('gobidata.csv', g.GeneExp(IDs,:)', '-append');
%             system('start ggobi gobidata.csv');
        end

        
        function Plot3(g, IDs, Group, Jitter, varargin);
            % PlotMatrix(IDs, Group, Jitter, varargin);
            %
            % makes a 3d plot of genes with specified IDs
            % IDs can be numbers or names but there must be 3 of them!
            % Jitter defaults to .3
            % varargin is passed to plot3
            
            if nargin<3 | isempty(Group)
                if ~isempty(g.Class)
                    Group = g.Class;
                else
                    Group = ones(g.nCells, 1);
                end
            end
            
            if g.nCells==0; error('Dude, no cells!'); end;

            [GroupNames, ~, Group] = unique(Group, 'stable');
            nGroups = max(Group);
            Colorset = HsvNotYellow(ceil(nGroups*1.2));
            
            if nargin<4
                Jitter = 0.3;
            end
            if length(IDs) ~= 3
                error('Plot3 needs 3 IDs!');
            end
                
            IDs = g.NamesToIDs(IDs);
            Sym = repmat('o+svp', [1 ceil(nGroups/5)]);

          
            cla; hold on;
            for i=1:nGroups
                MyCells = find(Group==i);
                nMy = length(MyCells);
                plot3(g.GeneExp(IDs(1),MyCells) + rand(1, nMy)*Jitter, ...
                    g.GeneExp(IDs(2),MyCells) + rand(1, nMy)*Jitter, ... 
                    g.GeneExp(IDs(3),MyCells) + rand(1, nMy)*Jitter, ...
                    Sym(i), 'color', Colorset(i,:), varargin{:});
            end
            hold off;
            
            xlabel(g.IDsToNames(IDs(1)), 'Interpreter', 'none');
            ylabel(g.IDsToNames(IDs(2)), 'Interpreter', 'none');
            zlabel(g.IDsToNames(IDs(3)), 'Interpreter', 'none');
            
            view([135 30]);
            
            grid on;
            legend(GroupNames);
        end

        function IDs = OrthogonalSet(g, n, StartSet, MaxLoops);
            % IDs = OrthogonalSet(n, MaxLoops)
            %
            % returns the IDs of n genes for which the sum of their cosine
            % angles is as small as possible.
            % 
            % StartSet gets you going (default is the two with lowest
            % angle)
            % MaxLoops (default 10) how many times to loop through looking
            % for stability before giving up
            
            if nargin<4
                MaxLoops = 10;
            end

            
            % compute matrix of cosine angles from L2 normalized vectors 
            fprintf('Computing kernel matrix ....');
            NormExp = bsxfun(@rdivide, g.GeneExp,sqrt(sum(g.GeneExp.^2,2)));
            C = NormExp*NormExp';
            fprintf('Done\n');
            
            if nargin<3 | isempty(StartSet)
                IDs = zeros(n,1);
                % find pair of minimum angle
                [dummy BestPair] = min(C(:));
                [IDs(1) IDs(2)] = ind2sub(g.nGenes*[1 1], BestPair);
                nPicked = 2;
                StartAt = 3;
            else
                IDs = zeros(n,1);
                IDs(1:length(StartSet)) = StartSet;
                StartAt = length(StartSet)+1;
                nPicked = length(StartSet);
            end
            
            StarterNames = g.IDsToNames(IDs(1:nPicked));
            fprintf('Starting with '); fprintf('%s ', StarterNames{:});
            fprintf('\n');


            for l=1:MaxLoops    
                AnyChanged = 0; % don't keep looping if nothing changed
                for i=StartAt:n
                    
                    fprintf('%d.%d: ', l, i);

                    OldID = IDs(i);
                    if l>1; fprintf('%s -> ', g.IDsToNames(OldID)); end;

                    % You can pick any gene that hasn't already been picked
                    OtherGenes = IDs(setdiff(1:nPicked,i));
                    Candidates = setdiff(1:g.nGenes,OtherGenes);
                    
                    % Compute sum of cosines with other picked cells
                    SumC = sum(C(OtherGenes,Candidates));
                    
                    % pick the gene with the smallest sum cosine angle
                    
                    [MinScore BestColumn] = min(SumC);
                    IDs(i) = Candidates(BestColumn);
                    
                    fprintf('%s has score %.5f\n', g.IDsToNames(IDs(i)), MinScore);
                    
                    nPicked = max(nPicked, i);
                    if IDs(i) ~= OldID; AnyChanged=1; end;
                end
                if AnyChanged==0; break; end;
                StartAt = 1;
            end
            
            AllWinners = g.IDsToNames(IDs);
            fprintf('Final genes: '); fprintf('%s ', AllWinners{:});
            fprintf('\n');

        end
        
        function [GeneNames, Scores] = Bisectors(g, Thresh, Reg, p, q, CandidateGenes);
            % [GeneNames, Scores] = Bisectors(Thresh, Reg, p, q, CandidateGenes);
            % search for a gene which, when it expresses >= Thresh, 
            % you find lots of other genes to express weakly. 
            % 
            % specifically for each candidate gene, it computes 
            % Ratio = (Mean1 + Reg)./(Mean0 + Reg) for all other genes
            % then sums 1./(1+(Ratio/q).^p) to get a score for each
            % candidate
            %
            % if Thresh is negative, it means this much times the mean 
            %
            % CandidateGenes is a list to try (default all)
            %
            % Output is an ordering of the Genes (best first), and their
            % ordered scores.
            
            if nargin<2; Thresh = -1; end
            if nargin<3; Reg = .1; end
            if nargin<4; p = 4; end
            if nargin<5; q = .1; end
            if nargin<6; 
                CandidateGenes = 1:g.nGenes;
            else
                CandidateGenes = g.NamesToIDs(CandidateGenes);
            end
            
            nc = g.nCells;
            ng = length(CandidateGenes);
            Stat = zeros(ng,1);
            nCells = zeros(ng,1);
            MyThresh = Thresh;
            for i=1:ng
                eg = CandidateGenes(i);
                if Thresh<0; 
                    MyThresh = -Thresh*mean(g.GeneExp(eg,:)); 
                end;
                Cells1 = find(g.GeneExp(eg,:)>=MyThresh);
                Cells0 = find(g.GeneExp(eg,:)<MyThresh);
                Mean1 = mean(g.GeneExp(:,Cells1),2);
                Mean0 = mean(g.GeneExp(:,Cells0),2);
                Ratio = (Mean1+Reg)./(Mean0+Reg);
                Logistic = 1./(1+(Ratio./q).^p);
                Stat(i) = sum(Logistic);
                nCells(i) = length(Cells1);
            end
                
                
%                 Mean1 = (sum(g.GeneExp(:,Cells),2)+Reg)/length(Cells);
%                 AllMean = (sum(g.GeneExp,2)+Reg)/nc;
%                 Markov = Mean1./AllMean;
%                 Stat(i) = sum(exp(-Markov/Thresh));
%                 Frac(i) = length(Cells)/nc;
%                 fprintf('%s: %d\n', g.IDsToNames(eg), Stat(i));
%             end
            
            Stat0 = Stat; Stat0(isnan(Stat))=0;
            [Scores GeneList] = sort(Stat0, 'descend');
            GeneNames = g.IDsToNames(CandidateGenes(GeneList));
            for i=1:min(10, ng)
                fprintf('%s has %f\n', GeneNames{i}, Scores(i));
            end
            
            clf; 
            semilogx(nCells, Stat, 'k.'); hold on; grid on
            h = text(nCells, Stat, g.GeneName(CandidateGenes), 'Interpreter', 'none');
            xlabel('n cells above threshold');
            ylabel('Score');
            set(h, 'color', 'b');

        end
            
                
        
%         function [Clu, GeneID, Thresh] = Bisect(g, CandidateGenes);
%             % [Clu, GeneID, Thresh] = Bisect(g, CandidateGenes);
%             %
%             % Finds a optimal bisection of the cells into two groups:
%             % GeneID<=Thresh and GeneID>Thresh
%             %
%             % Optimality criterion is low cosine angle of between group
%             % means
%             %
%             
%             % compute matrix of cosine angles from L2 normalized vectors 
%             % PREMATURE OPTIMIZATION IS THE ROOT OF ALL EVIL
% %             fprintf('Computing kernel matrix ....');
%              NormExp = bsxfun(@rdivide, g.GeneExp,sqrt(sum(g.GeneExp.^2,1)));
% %             C = NormExp*NormExp';
% %             fprintf('Done\n');
% 
%             % loop through candidate genes
%             CandidateGenes = g.NamesToIDs(CandidateGenes);
%             nCandidates = length(CandidateGenes);
%             Scores = zeros(nCandidates, g.nCells-1); % because there are nCells-1 possible nonempty splits
%             Sorteds = zeros(nCandidates, g.nCells);
%             for i=1:length(CandidateGenes)
%                 % sort matrix by expression of candidate gene
%                 [sorted order] = sort(g.GeneExp(CandidateGenes(i),:), 'ascend');
%                 Sorteds(i,:) = sorted;
%                 
%                 % cumulative sum of gene expression will help us make averages
%                 CumSum = cumsum(NormExp(:,order),2);
%                 
%                 % now mean of below threshold elements (leave out last to
%                 % avoid divide by zero)
%                 nBelow = 1:g.nCells-1;
%                 MeanBelow = bsxfun(@rdivide, CumSum(:,1:g.nCells-1), nBelow);
%                 % normalize
%                 nMeanBelow = bsxfun(@rdivide, MeanBelow, sqrt(sum(MeanBelow.^2,1)));
%                 
%                 % now mean of above threshold elemtns
%                 nAbove = g.nCells-1:-1:1;
%                 AboveSum = bsxfun(@minus, CumSum(:,end), CumSum);
%                 MeanAbove = bsxfun(@rdivide, AboveSum(:,1:g.nCells-1), nAbove);
%                 % normalize
%                 nMeanAbove = bsxfun(@rdivide, MeanAbove, sqrt(sum(MeanAbove.^2,1)));
%                 
%                 % now compute mean dot product of vectors with respective
%                 % means. Math trick lets us do this as |mu|^2*N
%                 
%                 Scores(i,:) = bsxfun(@times,sqrt(sum(MeanBelow.^2,1)),nBelow) + bsxfun(@times, sqrt(sum(MeanAbove.^2,1)), nAbove); 
%             end
%             
%             [BestScore, BestInd] = max(Scores(:));
%             [sBestGene, sBestThresh] = ind2sub([nCandidates, g.nCells-1], BestInd);
%             
%             GeneID = CandidateGenes(sBestGene);
%             Thresh = Sorteds(sBestGene,sBestThresh);
%             Clu = 1+(g.GeneExp(GeneID,:)>Thresh);
%             fprintf('Criterion: %s > %f, Score %f\n', g.IDsToNames(GeneID), Thresh, BestScore);
%         end
%                 
% 


                    
        
        function SortedNames = BestPredictors(g,nc, lambda, StartSet, MaxLoops)
            % BestPredictors(n, lambda, StartSet, MaxLoops)
            % returns the names of the n genes which best predict all others
            % in the data set via ridge regression with penalty lambda
            %
            % Startset gets you going with a list of user-specified genes
            % MaxLoops (default 1) means how many times to loop through
            % looking for stability, before giving up.
            
            if nargin<3; lambda = 1; end;
            if nargin<4; StartSet = {}; end;
            if nargin<5; MaxLoops = 10; end
            
            fprintf('Computing kernel matrix ....');
            X = g.GeneExp';
                
            if nargin<4 | isempty(StartSet)
                IDs = zeros(nc,1);
                Scores = zeros(nc,1);
                Kdata = X*X';
                Kpreds = zeros(g.nCells);
                %Kpreds = ones(g.nCells); % to deal with offsets we add a pseudogene of all ones
                StartAt = 1;
            else
                StartSet = g.NamesToIDs(StartSet);
                IDs = zeros(nc,1);
                Scores = zeros(nc,1);
                IDs(1:length(StartSet)) = StartSet;
                StartAt = length(StartSet)+1;
                
                Kpreds = X(:,StartSet)*X(:,StartSet)';
                Kdata = X*X' - Kpreds;
                %Kpreds = Kpreds + ones(g.nCells); % again, the pseudogene of all ones
            end
            
            
            
            fprintf(' Done\n');

            for l=1:MaxLoops    
                AnyChanged = 0; % don't keep looping if nothing changed
                for i=StartAt:nc
                    
                    fprintf('%d.%d: ', l, i);
                    % do we need to delete a gene?
                    if i<=nc & IDs(i)~=0
                        fprintf(' Losing %s: ', g.IDsToNames(IDs(i)));
                        v = X(:,IDs(i));
                        vv = v*v';
                        Kdata = Kdata+vv;
                        Kpreds = Kpreds-vv;
                    end
                    
                    % residual matrix
                    R = inv(eye(g.nCells) + Kpreds/lambda);
                    % loss
                    %L = trace(R*Kdata);
                    % total variance
                    Var = trace(Kdata);
                    % amount of unexplained variance
                    UnexpVar = trace(R*Kdata*R);
                    fprintf('prc var unexplained %f ', 100*UnexpVar/Var);

                    if i==nc+1
                        break;
                    end

                    DeltaL = zeros(g.nGenes,1);
                    %for gene=1:g.nGenes
                    parfor gene=1:g.nGenes
                        u = X(:,gene);
                        Ru = R*u;
                        uRu = u'*Ru;
                        % evaluate the effect of adding gene g
                        DeltaL(gene) = (-Ru'*Kdata*Ru + uRu^2)/(1+uRu);
                    end
                    [BestScore winner] = min(DeltaL);
                    WinnerName = g.IDsToNames(winner);
                    fprintf('Winner %s\n', WinnerName);

                    v = X(:,winner);
                    vv = v*v';
                    Kdata = Kdata-vv;
                    Kpreds = Kpreds+vv;
                    
                    if IDs(i)~= winner
                        AnyChanged = 1;
                    end
                    IDs(i)=winner;
                    Scores(i) = BestScore;

                end
                if AnyChanged==0; break; end;
            end
            
            %final score computation
            Kpreds = X(:,IDs)*X(:,IDs)';
%            Kpreds = X(:,IDs)*X(:,IDs)' + ones(g.nCells);
            Kdata = X*X' - Kpreds;
            R = inv(eye(g.nCells) + Kpreds/lambda);
            Var = trace(Kdata);
            UnexpVar = trace(R*Kdata*R);
            
            [SortedScores ScoreOrder] = sort(Scores, 'ascend');
            SortedIDs = IDs(ScoreOrder);
            
            SortedNames = g.IDsToNames(SortedIDs);
            if ~iscell(SortedNames); SortedNames = {SortedNames}; end; % to avoid crash with nPreds=1
            fprintf('Final genes: '); fprintf('%s ', SortedNames{:});
            fprintf('\nWith Scores: '); fprintf('%.1f ', SortedScores);
            fprintf('\nFinal prc var unexplained %f\n', 100*UnexpVar/Var);

        end
        
        function [Prediction, Weights] = Ridge(g, Target, Predictors, lambda)
            % [Prediction, Weights] = Ridge(g, Target, Predictors)
            % Target is gene to be predicted
            % Predictors is set of genes allowed to predict it
            %  - default is all but predicted gene
            % lambda is ridge penalty, default is 1000
            
            Target = g.NamesToIDs(Target);
            
            if nargin<3 | isempty(Predictors)
                Predictors = setdiff(1:g.nGenes, Target);
            else
                Predictors = g.NamesToIDs(Predictors);
            end
            
            if nargin<4
                lambda = 1000;
            end

            y = g.GeneExp(Target,:)';
            X = g.GeneExp(Predictors,:)';
            %X = [g.GeneExp(Predictors,:)', ones(g.nCells, 1)];
            K = X*X';
            
            % compute prediction using hat matrix
            H = eye(g.nCells) - lambda*inv(lambda*eye(g.nCells) + K);
            Prediction = H*y;
            
            VarExp = 1 - var(Prediction-y)/var(y);
            
            % now compute weights
            Weights = (y-Prediction)'*X/lambda;
            
            subplot(1,2,1); cla; hold on
            g.Scatter(Prediction, Target);
            plot(xlim, xlim, 'k');
            xlabel('Prediction');
            title(sprintf('%f percent variance explained', 100*VarExp));
            
            subplot(1,2,2);
            barh(Weights); 
            set(gca, 'yticklabel', g.IDsToNames(Predictors));
            xlabel('Weight');
            
%             [b,se,pval,inmodel] = stepwisefit(g.GeneExp(Predictors,:)', g.GeneExp(Target,:)');
%             Weights = b.*inmodel(:);
%             Prediction = g.GeneExp(Predictors,:)'*Weights;
        end
        
        function gOut = Residual(g, Predictors, lambda)
            % gOut = Residual(g, Predictors, lmabda)
            %
            % Takes out a common effect of predictor genes on the gene set
            % by fitting a linear regression (no offset) then subtracting
            % this prediction
            %
            % lambda is ridge term - make it high!
            %
            % Use this to remove contamination by having predictors you
            % know are not really expressed in your cells

            Predictors = g.NamesToIDs(Predictors);
            
            % transpose matrices - because this is much easier to understand
            Y = g.GeneExp';
            X = g.GeneExp(Predictors,:)';
            
            % do ridge regression with trick of adding fake observations
            Yplus = [Y; zeros(length(Predictors), g.nGenes)];
            Xplus = [X; lambda*eye(length(Predictors))];


            % find W s.t Y~X*W
            W = X\Y;
            Yhat = X*W;
            
            gOut = g;
            gOut.GeneExp = g.GeneExp - (X*W)';
            
        end


        function Lasso(g, Target);
            % LassoGene(g, Target);
            % predict a target by lasso from all other genes in the set
            %

            Target = g.NamesToIDs(Target);

            y = g.GeneExp(Target,:)';
            Predictors = setdiff(1:g.nGenes, Target);
            X = g.GeneExp(Predictors ,:)';

            [b stats] = lasso(X, y, 'NumLambda', 100, 'PredictorNames', g.GeneName(Predictors), 'CV', 2, 'DFMax', 20);

            lassoPlot(b, stats, 'PredictorNames', g.GeneName(Predictors), 'PlotType', 'Lambda');
            lassoPlot(b, stats, 'PredictorNames', g.GeneName(Predictors), 'PlotType', 'CV');
        end
        
        function v = On(g, Gene, Thresh)
            % returns a binary vector containing 1 for cells when Gene's
            % expression is at least Thresh (default 5)
            % 
            % If Thresh is negative, threshold is -Thresh* mean expression
            
            Gene = g.NamesToIDs(Gene);
            
            if nargin<3
                Thresh = 5;
            end
            if Thresh<0
                Thresh = -Thresh*mean(g.GeneExp(Gene,:));
            end
            
            v = g.GeneExp(Gene,:)>=Thresh;
        end
        
        function Selected = Friends(g, Gene, nFriends, Group, Rank);
            % Selected = Friends(Gene, nFriends, Group, Rank);
            % 
            % find genes best correlated with target Gene: can be name, ID,
            % or numerical vector of length g.nCells. (Actually shows
            % enemies first, unless nFriends is negative)
            %
            % if you don't give specify an output, it will interactively
            % make scatter plots of the genes most positively and negatively 
            % correlated with a gene of interest. Press return to see the
            % next, press a 'w' then return to launch Allen and OMIMM info 
            % in a web browser
            %
            % Optional input Rank (default 2) ranks the data first (if 1 doesn't adjust
            % for ties, to save time; if 2 deals with ties).
            
            if nargin<3 | isempty(nFriends), nFriends = 10; end;
            if nargin<4 | isempty(Group)
                if ~isempty(g.Class)
                    Group = g.Class;
                else
                    Group = []; 
                end
            end
            if length(Group)==1;
                Group = [];
            end
            if nargin<5
                Rank = 2;
            end
            
            if g.nCells<1
                warning('Dude, no cells.');
                return
            end
            
            
            
            % get ranks. Because tiedranks takes a long time, option 1 just
            % uses sorting.
            if Rank==1
                SmallNumber = 1e-3;
                Data = RankOrder(g.GeneExp', SmallNumber);
                Reg = 2;
            elseif Rank==2 % if you have more time, use option 2
                Data = tiedrank(g.GeneExp');
                Reg=5;
            else
                Data = g.GeneExp';
                Reg = 2;

            end

            if length(Gene)==g.nCells % did you pass a numerical vector for target?
                StartFrom = 1;
                if Rank==1
                    Target = RankOrder(Gene(:), SmallNumber);
                elseif Rank==2
                    Target = tiedrank(Gene(:));
                else
                    Target = Gene(:);
                end
            elseif ismember(Gene, g.GeneName)
                StartFrom = 1;
                Target = Data(:,g.NamesToIDs(Gene));
            elseif ismember(Gene, fieldnames(g.CellInfo))
                StartFrom = 1;
                Target = getfield(g.CellInfo,Gene);
            else
                warning(sprintf('Could not find Gene or Cell ID %s', Gene));
                return;
            end
            
            % compute correlation coefficients with BaseGene
            % regularize the z-scores
            %z = zscore(Data);
            z = bsxfun(@rdivide,bsxfun(@minus, Data, mean(Data)), std(Data)+Reg);
            %zTarget = zscore(Target(:)); 
            zTarget = (Target(:)-mean(Target(:)))/std(Target(:)+Reg); 

            c = zTarget'*z/double(g.nCells);

            [dsorted dorder] = sort(c, 'descend');
            [asorted aorder] = sort(c, 'ascend');
            if nFriends>0
                sorted = [asorted(1:nFriends), dsorted(StartFrom:nFriends+1)]; % first dsorted will be gene itself
                order = [aorder(1:nFriends), dorder(StartFrom:nFriends+1)];
            else
                % nFriends negative means show positive correlations first
                sorted = [dsorted(StartFrom:-nFriends+1), asorted(1:-nFriends)]; % first dsorted will be gene itself
                order = [dorder(StartFrom:-nFriends+1), aorder(1:-nFriends)];
            end
            

            Selected = [];
            for i=1:length(order)
%                 figure(1); clf
                g.Scatter(Gene, order(i), Group);
%                 
                set(gca, 'FontSize', 12);
                drawnow;

                %web(['https://en.wikipedia.org/wiki/Special:Search/' GeneName{order(i)}]);

                fprintf('%s: Correlation %f.', g.IDsToNames(order(i)), sorted(i));
                str = input('','s');
                
                if strcmp(str, 'q')
                    break;
                end
                if ~isempty(str)
                    Selected = [Selected, order(i)];
                end

                if strcmp(str, 'w')
                    system(['start chrome "http://mouse.brain-map.org/search/show?search_type=gene&search_term=' g.GeneName{order(i)} '"']);
                    system(['start chrome "http://www.ncbi.nlm.nih.gov/omim/?term=' g.GeneName{order(i)} '"']);
                end

            end
            
            if ~isempty(Selected)
                Selected = g.IDsToNames(Selected);
            elseif nargout>0
                Selected = g.IDsToNames(order);
                
            end

        end
        
        function BoxPlot(g, Gene, Class)
            % BoxPlot(Gene, Class)
            % does a box plot to show how much a gene is expressed in each
            % class. If Gene is an array, it does all of them.
            %
            % if Class is missing, uses g.Class
            %
            % if Gene is a vector of length g.nCells, plots these values
            
            if nargin<3
                Class = g.Class;
            end
            
            if length(Gene)==g.nCells
                boxplot(Gene, Class);
                title('User data');
                set(gca, 'xticklabelrotation', 45);
                ylim([0 max(ylim)]);
            else
                if isstr(Gene) & Gene(end)=='*'
                    Gene = g.NameStartsWith(Gene(1:end-1));
                else
                    Gene = g.NamesToIDs(Gene);
                end
                for i=1:length(Gene(:)')
                    clf; hold on
                    boxplot(g.Exp(Gene(i)), Class);
                    ColNames = get(gca, 'xticklabel');
                    for j=1:length(ColNames)
                        mu = mean(g.CellSubset(ColNames{j}).Exp(Gene(i)));
                        plot(j, mu, 'kx');
                    end
                    title(g.IDsToNames(Gene(i)));
                    set(gca, 'xticklabelrotation', 45);
                    ylim([0 max(ylim)]);
                    if i<length(Gene)
                        str = input('', 's');
                        if strcmp(str, 'w')
                            system(['start chrome "http://mouse.brain-map.org/search/show?search_type=gene&search_term=' g.GeneName{Gene(i)} '"']);
                            system(['start chrome "http://www.ncbi.nlm.nih.gov/omim/?term=' g.GeneName{Gene(i)} '"']);
                        end
                    end
                end
            end
            
        end
        
        function Heatmap(g, Genes, Class, ClassNames)
            % Heatmap(g, Genes, Classes)
            %
            % Plot a heatmap of a subset of genes showing the mean
            % expression of each class
            %
            % Genes: subset of Genes to show
            % Class: class ID for each gene
            
            if nargin<3; Class = g.Class; end;
            Genes = g.IDsToNames(g.NamesToIDs(Genes)); % cut any missing ones
            
            [c, ClassNames] = grp2idx(Class);
            nClasses = length(ClassNames);
            nGenes = length(Genes);
            ColorMat = zeros(nClasses, nGenes);

            % make expression matrix
            Exp = g.Exp(Genes);
            for i = 1:nClasses
                ColorMat(i,:) = mean(log(1+Exp(:,c==i)),2);
            end
            
            % sort genes into the best order
            [~,BestClass] = max(ColorMat,[], 1);
            [~,GeneOrder] = sort(BestClass, 'ascend');

            clf;

            % draw colormap and label Classes
            imagesc(ColorMat(:,GeneOrder)');
            set(gca, 'ytick', [])
            set(gca, 'xtick', 1:nClasses);
            set(gca, 'xticklabel', ClassNames);
            set(gca, 'xticklabelrotation', 90);
            set(gca, 'position', [0.18 0.145, 0.63 0.8]);
            hold on

            % now draw gene names, putting them to the side to squeeze in loads
            % plot shape parameters:
            l = .5;
            r = nClasses+.5;
            s1 = nClasses/50;
            s2 = nClasses/10;
            fsz = 9;
            gp = 0.05;
            colormap(hot);
            %caxis([-1 1]);

            for i=1:4:nGenes
                text(l-s1-gp, i, Genes{GeneOrder(i)}, 'HorizontalAlignment', 'Right', 'fontsize', fsz);
                h=plot([l-s1 l+s1], [i i], 'k');
                set(h, 'clipping', 'off');
            end

            for i=2:4:nGenes
                text(r+s1+gp, i, Genes{GeneOrder(i)}, 'HorizontalAlignment', 'Left', 'fontsize', fsz);
                h=plot([r-s1 r+s1], [i i], 'k');
                set(h, 'clipping', 'off');
            end

            for i=3:4:nGenes
                text(l-s2-gp, i, Genes{GeneOrder(i)}, 'HorizontalAlignment', 'Right', 'fontsize', fsz);
                h=plot([l-s2 l+s1], [i i], 'k');
                set(h, 'clipping', 'off');
            end

            for i=4:4:nGenes
                text(r+s2+gp, i, Genes{GeneOrder(i)}, 'HorizontalAlignment', 'Left', 'fontsize', fsz);
                h=plot([r-s1 r+s2], [i i], 'k');
                set(h, 'clipping', 'off');

            end

        end

        
        function [gFactors, gScores] = Nnmf(g, SeedGenes)
            % [gFactors gScores] = Nnmf(g, m);
            % non-negative matrix factorization with m factors
            %
            % so g.GeneExp ~ gFactors.GeneExp * gScores.GeneExp
            %
            % thus gFactors has nGenes genes and m pseudocells
            % gScores has m pseudogenes and nCells cells

%           [m1 m2] = nnmf(g.GeneExp, k);

            SeedGenes = g.NamesToIDs(SeedGenes);
            k = length(SeedGenes);
                        
            % compute g.GeneExp = m1 (nGenes x k) * m2 (k x nCells)
            MaxIter = 100;
            G = g.GeneExp;
            h0 = G(SeedGenes,:);
            w0 = max(0, G/h0);
            
            tol=1e-8;
            Err = inf;
            for i=1:MaxIter
                OldErr = Err;
                
                if 0
                    w = max(0, G/h0);
                    h = max(0, w\G);
                else                
                    numer = w0'*G;
                    h = h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer)));
                    numer = G*h';
                    w = w0 .* (numer ./ (w0*(h*h') + eps(numer)));
                end

                
                Resid = G-w*h;
                Err = sum(Resid(:).^2);
                fprintf('Iter %d, error %f\n', i, Err);
                if (OldErr-Err)/OldErr<tol; break; end
                w0 = w;
                h0 = h;
            end
            
            gFactors = GeneSet;
            gFactors.GeneExp = w0;
            gFactors.GeneName = g.GeneName;
            gFactors.nGenes = g.nGenes;
            gFactors.nCells = k;
            
            gScores = GeneSet;
            gScores.GeneExp = h0;
            for i=1:k
                gScores.GeneName{i} = sprintf('Factor %d', i);
            end
            gScores.nGenes = k;
            gScores.nCells = g.nCells;
            gScores.CellInfo = g.CellInfo;
        end
        
        function [gFactors, gScores] = PCA(g, m)
            % [gFactors, gScores] = g.PCA(m)
            % PCA with m factors (subtracting mean expression of each gene)
            %
            % so g.GeneExp ~ gFactors.GeneExp * gScores.GeneExp
            %
            % thus gFactors has nGenes genes and m pseudocells
            % gScores has m pseudogenes and nCells cells
            %
            % You might want to run Log first
            
            % subtract mean from each gene
            MeanExp = mean(g.GeneExp,2);
            Exp0 = bsxfun(@minus, g.GeneExp, MeanExp);
            
            % to compute PCs, compute eigs of cells covariance (much quicker)
            Cov = Exp0' * Exp0;
            [v, d] = eigs(double(Cov), m); % each column of v is a score, normalized
            Scores = sqrt(d)*v';
            Factors = Exp0/Scores;
            
            gFactors = GeneSet;
            gFactors.GeneExp = Factors;
            gFactors.GeneName = g.GeneName;
            gFactors.nGenes = g.nGenes;
            gFactors.nCells = m;
            
            gScores = GeneSet;
            gScores.GeneExp = Scores;
            gScores.nGenes = m;
            gScores.nCells = g.nCells;
            gScores.CellInfo = g.CellInfo;
            gScores.CellName = g.CellName;
            gScores.Class = g.Class;

            for i=1:m
                gScores.GeneName{i} = sprintf('Factor %d', i);
                gFactors.CellName{i} = sprintf('Factor %d', i);
            end

        end
        
        function [gFactors, gScores] = NBpca(g, m, r)
            % [gFactors, gScores] = g.NBpca(m)
            % 
            % Does negative-binomial PCA with m factors 
            % (not including constant, default m=1)
            %
            % so g.GeneExp~NB(exp(gFactors.GeneExp * gScores.GeneExp),r)
            % (Default r=2) 
            %
            % output gFactors has nGenes genes and m+1 pseudocells
            % gScores has m+1 pseudogenes and nCells cells
            % 
            % if no outputs, it plots factor 1 vs factor 2 (which could be
            % offset)
            % 
            % note that you should first scale then take a subset of 
            % high-expressing genes e.g. g.ScaleCell(1).GoodExpressers(5,2).NBpca;
            
            if nargin<2; m=1; end
            if nargin<3; r=2; end
            
            [w, x] = NBpca(double(g.GeneExp'), r, m, [], [], .1, .1);
            % w is nCells by m+1
            % x is m+1 by nGenes
            
            if nargout==0
                clf; hold on
                yCoord = w(2,:);
                plot(w(1,:), yCoord, 'b.');
                hand = text(w(1,:), yCoord, g.GeneName); 
                set(hand, 'color', 'b');
                grid on
                drawnow
            end
            
            gFactors = GeneSet;
            gFactors.GeneExp = w';
            gFactors.GeneName = g.GeneName;
            gFactors.nGenes = g.nGenes;
            gFactors.nCells = m+1;
            
            gScores = GeneSet;
            gScores.GeneExp = [x'; ones(1,g.nCells)];
            gScores.nGenes = m+1;
            gScores.nCells = g.nCells;
            gScores.CellInfo = g.CellInfo;
            gScores.CellName = g.CellName;
            gScores.Class = g.Class;
            
            for i=1:m
                gScores.GeneName{i} = sprintf('Factor %d', i);
                gFactors.CellName{i} = sprintf('Factor %d', i);
            end
            gScores.GeneName{m+1} = 'Constant';
            gFactors.CellName{m+1} = 'Offset';


        end
            
            
            
        function [gFactors, gScores] = FactorAn(g, m)
            % [gFactors gScores] = g.FactorAn(m);
            % factor analysis with m factors
            %
            % so g.GeneExp ~ gFactors.GeneExp * gScores.GeneExp
            %
            % thus gFactors has nGenes genes and m pseudocells
            % gScores has m pseudogenes and nCells cells
            
            [lambda psi T stats F] = factoran(g.GeneExp', m);
            
            gFactors = GeneSet;
            gFactors.GeneExp = lambda;
            gFactors.GeneName = g.GeneName;
            gFactors.nGenes = g.nGenes;
            gFactors.nCells = m;
            
            gScores = GeneSet;
            gScores.GeneExp = F';
            for i=1:m
                gScores.GeneName{i} = sprintf('Factor %d', i);
            end
            gScores.nGenes = m;
            gScores.nCells = g.nCells;
            gScores.CellInfo = g.CellInfo;
        end
        
        function gFactors = PredictFrom(g, gScores)
            % given a subset of genes (or pseudogenes), predict the rest
            %
            % That is, given a gene set g, and a m x nCells set of pseudogenes 
            % find gFactors with nGenes and m pseudocells such that 
            % 
            % g.GeneExp ~ gFactors.GeneExp * gScores.GeneExp
            
            gFactors = GeneSet;
            gFactors.GeneExp = g.GeneExp / gScores.GeneExp;
            gFactors.GeneName = g.GeneName; 
            gFactors.nGenes = g.nGenes;
            gFactors.nCells = gScores.nGenes;
            gFactors.CellName = gScores.GeneName;
        end
        
        function gCorrels = Correlations(g, g2)
            % given two sets of genes (or pseudogenes), find the correlation
            % coefficients of each pair of genes
            %
            % That is, given a m x nCells gene set g, and a m2 x nCells
            % gene set g2, find the m x m2 matrix of correlation
            % coefficients
            
            z = zscore(g.GeneExp,1,2);
            z2 = zscore(g2.GeneExp,1,2);
            c = z*z2'/g.nCells;

            gCorrels = GeneSet;
            gCorrels.GeneExp = c;
            gCorrels.GeneName = g.GeneName; 
            gCorrels.nGenes = g.nGenes;
            gCorrels.nCells = g2.nGenes;
            gCorrels.CellName = g2.GeneName;
        end



        
        function Biplot(gFactors, gScores, ScoreColors, FactorColors, varargin)
            % gFactors.Biplot(gScores, ScoreColors, FactorColors, varargin)
            % 
            % do a biplot for factor analysis/nnmf/etc
            %
            % it will plot arrows and gene names for gFactors
            % optional input gScores will make it show points also (use [] not to show them)
            % 
            % FactorColors and ScoreColors will color arrows and cells
            %
            % varargin is passed onto biplot (e.g. 'positive', true);
            
            if nargin<2
                gScores = [];
            end
            
            MaxFactors = 500; % don't plot more than this - you can't read them anyway
            TotFact = sum(zscore(gFactors.GeneExp).^2,2);
            [~, order] = sort(TotFact, 'descend');
            gFactors = gFactors.GeneSubset(order(1:min(MaxFactors,gFactors.nGenes)));
            
            if ~isempty(gScores)
                handle = biplot(gFactors.GeneExp, 'scores', gScores.GeneExp', ...
                    'varlabels', gFactors.GeneName, varargin{:});
            else
                handle = biplot(gFactors.GeneExp, 'varlabels', gFactors.GeneName, varargin{:});
            end
            
            % axis labels if possible
            if ~isempty(gFactors.CellName)
                xlabel(gFactors.CellName{1});
                ylabel(gFactors.CellName{2});
                if gFactors.nCells>2
                    zlabel(gFactors.CellName{3});
                end
            end
            
            uistack(handle(1:gFactors.nGenes), 'bottom');
            
            % colors for scores
            if nargin<3 & ~isempty(gScores);
                ScoreColors = gScores.Class;
            else 
                ScoreColors = [];
            end
            if ~isempty(ScoreColors) & ~isempty(gScores)
                [c ia ic] = unique(ScoreColors);
                Colorset = HsvNotYellow(ceil(max(ic)*1.2));
                Color = Colorset(ic,:);
                for i=1:length(ic)
                    set(handle(gFactors.nGenes*3 + i), 'color', Color(i,:));
                end
            end
            
            % colors for factors
            if nargin>=4 & ~isempty(FactorColors) 
                [c ia ic] = unique(FactorColors);
                Colorset = hsv(max(ic)+1);
                Color = Colorset(ic,:);
                for i=1:length(ic)
                    set(handle(i + [0 1 2]*gFactors.nGenes), 'color', Color(i,:));
                end
            end
            
            if gFactors.nCells>2
                view([135 45]);
            end
        end
        
        function c = Cosine(g, Cells);
            % return a cosine similarity for all cells with target Cell
            % if Cells is a vector of length nGenes, it uses that as target.
            % if it is a vector of a different length, it averages those
            % cells
            
            if length(Cells) == g.nGenes
                Mean = Cells; 
            else
%                 Cells = g.NamesToIDs(Cells);
                Mean = mean(g.GeneExp(:,Cells),2);
            end
            Mean = Mean/norm(Mean);
            
            % norm of each cell in the database, to take dot prod
            Norm = sqrt(sum(g.GeneExp.^2, 1));
            
            c = Mean' * bsxfun(@rdivide,g.GeneExp,Norm);
        end
            
        
        function Neighbors = NearestNeighbors(g, Cells, nNeighbors);
            % Neighbors = NearestNeighbors(Cells, nNeighbors)
            %
            % sorts cells in order of cosine angle distance to the mean of
            % the Cells provided, and returns an array with 1 if each cell
            % is in the nearest nNeighbors, 0 otherwise.
            % if Cells is a vector of length nGenes, it uses that as mean.
            
           
%             if nargin<3
%                 nNeighbors = g.nCells;
%             end
            
            if length(Cells) == g.nGenes
                Mean = Cells; 
            else
                Mean = mean(g.GeneExp(:,Cells),2);
            end
            
            % norm of each cell in the database, to take dot prod
            Norm = sqrt(sum(g.GeneExp.^2, 1));
            
            Cosine = Mean' * bsxfun(@rdivide,g.GeneExp,Norm);
            [sorted order] = sort(Cosine, 'descend');
            Neighbors = (Cosine>=sorted(nNeighbors));
        end
        
        function CellSimilarityPlot(g, c1, c2);
            % CellSimilarityPlot(c1, c2);
            % plots something to show whicch genes contribute to similarity
            % via cosine angle: 
            %
            % c1 and c2 can be cell names, IDs, or vectors of length nGenes
            %
            % 
            
            if isstr(c1)
                c1 = find(strcmp(c1, g.CellName));
            end
            
            if isstr(c2)
                c2 = find(strcmp(c2, g.CellName));
            end
            
            if length(c1)==1 
                v1 = g.GeneExp(:,c1);
            else
                v1 = c1;
            end
            v1 = v1./norm(v1);
            
            if length(c2)==1
                v2 = g.GeneExp(:,c2);
            else
                v2 = c2;
            end
            v2 = v2./norm(v2);
            
            clf; subplot(1,2,1); hold on;
            % draw contours of constant product
            xr = .001:.001:1;
            m = nanmax([v1(:); v2(:)]);
            for r=m^2.*2.^(-10:0);
                loglog(xr, r./xr, 'k:');
            end
            
            loglog(v1,v2, 'b.');
            %set(gca, 'xscale', 'log'); set(gca, 'yscale', 'log');
            
            h = text(v1*1.01, v2, g.GeneName, 'Interpreter', 'none');
            
            for i=1:length(h)
                set(h(i), 'ButtonDownFcn', {@MDSCallback, g.GeneName{i}});
            end
            set(h, 'color', 'b');
            
            axis([0 m 0 m]);
            
            subplot(1,2,2);
            semilogx(1:g.nGenes, cumsum(sort(v1.*v2, 'descend')));
            xlabel('number of genes');
            ylabel('cumulative dot product');
            grid on
            
            fprintf('Cosine angle %f\n', v1(:)'*v2(:));
            
        end
        
        function CellDifferencePlot(g, c1, c2, Reg);
            
            if isstr(c1)
                c1 = find(strcmp(c1, g.CellName));
            end
            
            if isstr(c2)
                c2 = find(strcmp(c2, g.CellName));
            end

            if nargin<4
                Reg = 1;
            end
            
            if length(c1)==1
                v1 = g.GeneExp(:,c1);
            else
                v1 = c1;
            end
            
            if length(c2)==1
                v2 = g.GeneExp(:,c2);
            else
                v2 = c2;
            end

            clf; 
            x = (v1+v2)/2;
            y = (v1+Reg)./(v2+Reg);
            loglog(x, y, '.');
            h = text(x*1.01, y, g.GeneName, 'Interpreter', 'none');
            grid on
            set(h, 'color', 'b');
            xlabel('Mean Expression');
            ylabel('Expression Ratio (Gp1 / Gp2)');
        end

        
        function h = GeneScatter(g, c1, c2, GenesToShow, Jitter, Highlight)
            % h = GeneScatter(g, c1, c2, GenesToShow)
            % plots expression of all genes in c1 vs in c2
            % with text labels by points.            
            %
            % If c1 and c2 are vectors of length g.nGenes, it plots these
            % Otherwise then can be groups of cells: it uses IdentifyCells, and
            % takes the 34% trimmean of those cells. c2 empty or missing => all but
            % c1. If some cells are in both, they are removed from group c2
            %
            % GenesToShow says which genes to plot. If this is a single number, 
            % the function will plot that many that have largest difference index.
            % If it is an array, plot those genes. Default: top 100.
            % 
            % returns vector of handles
            % 
            % Jitter is a number, how much to jitter (default 0.02)
            % Highlight is a list of genes to plot in red
            
            reg = 1; % regularization parameter, for choosing which to show
            
            % names will be useful to plot stuff later
            oc1 = c1; 
            if nargin>=3; 
                oc2 = c2; 
            else
                oc2 = '';
            end;
            
            % compute what genes to show - default 100 top of each class
            if nargin<4; GenesToShow = min(100,g.nGenes); end;
            if nargin<5; Jitter = .02; end;
            if nargin<6; 
                Highlight = [];
                RegularColor = 'b';
                HighlightColor = 'k';
            else
                RegularColor = [.4 .4 1];
                HighlightColor = 'k';
            end; 
            Highlight = g.NamesToIDs(Highlight);

            
            % get mean vector for c1
            if isnumeric(c1) & (length(c1)==g.nGenes);
                % numeric vector goes right through
                v1all = c1;
            else
                % identify from groups
                c1 = g.IdentifyCells(c1);
                v1all = trimmean(g.GeneExp(:,c1),34, '', 2);
                
                % if no second argument provided, it is the rest
                if nargin<3 | isempty(c2)
                    c2 = setdiff(1:g.nCells,c1);
                end
                
            end
            
            
            % get vector for c2
            if isnumeric(c2) & length(c2) == g.nGenes;
                % numeric vector goes right through
                v2all = c2;
            else
                % identify from groups and remove any from c1
                c2 = setdiff(g.IdentifyCells(c2), c1);
                v2all = trimmean(g.GeneExp(:,c2),34, '', 2);
            end
            
            
            if length(GenesToShow)==1
                d = (v1all-v2all)./(v1all+v2all+reg);
                [~, order1] = sort(abs(d), 'descend');
                [~, order2] = sort(abs(d), 'ascend');
                GenesToShow = [order1(1:GenesToShow); order2(1:GenesToShow)];
            end

            v1 = v1all(GenesToShow);
            v2 = v2all(GenesToShow);


            % compute mean expression of these cell sets
            
            if length(v1)>5000
                s = sprintf('Do you really want to plot %d genes?', length(v1));
                resp = input(s, 's');
                if isempty(resp) | ~ismember(resp(1), 'yY'); return; end;
            end
            
            GeneNames = g.IDsToNames(GenesToShow);
            
            cla; hold on
            x = v1 + Jitter*rand(size(v1));
            y = v2 + Jitter*rand(size(v2));
            BlueGenes = ~ismember(GenesToShow,Highlight);
            h1 = text(double(x(BlueGenes))*1.01, double(y(BlueGenes)), GeneNames(BlueGenes), 'Interpreter', 'none');
            set(h1, 'color', RegularColor);
            plot(x(BlueGenes), y(BlueGenes), '.', 'color', RegularColor);
            
            RedGenes = ismember(GenesToShow,Highlight);
            plot(x(RedGenes), y(RedGenes), '.', 'color', HighlightColor);
            h2 = text(double(x(RedGenes))*1.01, double(y(RedGenes)), GeneNames(RedGenes), 'Interpreter', 'none');
            grid on
            set(h2, 'color', HighlightColor);
            h = [h1; h2];
            
            MaxVal = max([v1(:); v2(:)]);
            axis([0 1 0 1]*MaxVal);
            plot([0 MaxVal], [0 MaxVal], 'k:');
            
            if isstr(oc1); xlabel(oc1, 'Interpreter', 'none'); else; xlabel('user'); end
            if isstr(oc2); ylabel(oc2, 'Interpreter', 'none'); else; ylabel('user'); end
            
            for i=1:length(h)
                set(h(i), 'ButtonDownFcn', {@MDSCallback, GeneNames{i}});
            end

            % avoid output spew if no output
            if nargout==0; clear h; end
        end
            
        
        function [g, Decisions] = DivisiveClustering(g, SplitClass, Params, Recurse)
            % g = DivisiveClustering(SplitClass, Params) 
            % Recursive program to do divise clustering of cells
            %
            % Those cells whose entry in g.Class exactly match SplitClass will be
            % split into two groups by calling Params.SplitFn(). Then g.Class will
            % be updated for these cells, appending ".Gene1" and ".Gene2"
            % where Gene1 and Gene2 stand for the genes of highest choice
            % probability for the two subsets.
            %
            % After this, DivisiveClustering will be called recursively for
            % these two subclasses.
            %
            % if you call with SplitClass empty or missing, it will wipe
            % g.Class. So be careful!
            %
            % Params is a structure containing a field SplitFn saying which 
            % split function to call, and other parameters that get passed
            % to it. Use DefaultSplitParams to create default
            %   
            %
            % returns g with Class added
            % Second optional output Decisions is a structure array with
            % one entry per split, containing:
            %   w: 2 by nGenes array of weights with NaN for genes not considered
            %   y: 2 by nCells array of class membership scores with NaN for cells not considered
            %   ParentClass: string saying who was split
            %   ChildClass1, ChildClass2: strings saying what it got split into
            
%             
%             if nargin<3
%                 DisplayGenes={};
%             else
%                 DisplayGenes = intersect(DisplayGenes, g.GeneName);
%             end

            if nargin<3 
                Params = DefaultSplitParams;
            end
            
            if nargin<4
                Recurse=0;
            end
            
            
            g.GeneExp = double(g.GeneExp); % because MATLAB
            
            % is this the top-level call? if so, wipe existing 
            if Recurse==0
                SplitClass = '';
                MyCells = 1:g.nCells;
                g.Class = cell(g.nCells,1);
                g.Class(1:g.nCells) = cellstr(SplitClass);
            else
                MyCells = find(strcmp(SplitClass, g.Class));
            end
            
            fprintf('Splitting %s\n', SplitClass);
            
            

            [Clu, Suffix1, Suffix2, Decision] = Params.SplitFn(g,MyCells, Params, SplitClass,0);
            %[Clu, Suffix1, Suffix2, Decision] = g.Split(MyCells, Params, SplitClass,0);
            %[Clu, Suffix1, Suffix2, Decision] = SplitLP(g,MyCells, Params, SplitClass,0);
            Decision.ParentClass = SplitClass;
               
            Name1 = [SplitClass, '.', Suffix1];
            Name2 = [SplitClass, '.', Suffix2];
            Decision.ChildClass1 = Name1;
            Decision.ChildClass2 = Name2;

            if ~isempty(Clu)
                % rename split classes
                
                c1 = find(Clu==1);
                c2 = find(Clu==2);
                g.Class(MyCells(c1)) = cellstr(Name1);
                g.Class(MyCells(c2)) = cellstr(Name2);
                
                % now recusively split each one again
                if Decision.UserResponse~='f'
                    [g, Decision1] = g.DivisiveClustering(Name1, Params, 1);
                    [g, Decision2] = g.DivisiveClustering(Name2, Params, 1);
                else
                    Decision1 = [];
                    Decision2 = [];
                end

                % find genes most strongly expressed in each half
%                 ExpRat = (mean(g1.GeneExp,2)+eps)./(mean(g2.GeneExp,2)+eps);
%                 [dummy, TopFor1] = max(ExpRat);
%                 [dummy, TopFor2] = min(ExpRat);
%                 
%                 Suffix1 = g.IDsToNames(TopFor1);
%                 Suffix2 = g.IDsToNames(TopFor2);
                Decisions = [Decision, Decision1, Decision2];
            else
                Decisions = Decision;
            end
            
            if Recurse==0
                % sort by class name, only if top iteration
                [~, order] = sort(g.Class);
                g = g.CellSubset(order);
            end
            
        end
        
        function [Clu, Suffix1, Suffix2, Decision, BestScore] = Split(gAll, TheseCells, DisplayGenes, PlotTitle, ScoreOnly)
            % split into an L shape
            % output 1 or 2 for each cluster, 0 for neither
            % output [] not to split
            % optional outputs Suffix1 and Suffix2 are for saying what to
            % call them (used in divisive hierarchical clustering).
            % Decision is a structure with fields:
            %   w: 2 x gAll.nGenes array of weights
            %   y: gAll.nCells array of scores
            %   they will be NaN for genes or cells that were not
            %   considered.
            %
            % FinalScore is the score of the chosen split ()
            
            if nargin<2 | isempty(TheseCells)
                g = gAll;
            else
                g = gAll.CellSubset(TheseCells);
            end
            
            if nargin<3; DisplayGenes = []; end
            if nargin<4; PlotTitle = []; end
            if nargin<5; ScoreOnly = 1; end

            
            % first find VectorLength most bimodal genes
            VectorLength = 100;
            nStarts = VectorLength;
            alpha = 0.05;
            %alpha = 10;
            
            beta = 1/g.nCells;
            p=2;
            d=2;
            b=1/d;
            Jitter = .25;
            Thresh = ((b+1)/2).^p; % corresponds to x=1, y=0
            Thresh = -inf; % CLASSIFY EVERYTHING!!!
            
            SkipThresh = 1e-3;

            nTops = 60; % how many to show in colored words

            
            
            pow = 2;


            
            fprintf('Finding %d most bimodal...', VectorLength);
            Bim = g.Bimodality;
            [BimSorted, BimOrder] = sort(Bim, 'descend');
            gBim = g.GeneSubset(BimOrder(1:VectorLength));
            fprintf('Done\n');
            
            
            %now do divisive clustering on them
            X = gBim.ScaleGene(1).GeneExp;

            % parameters
            
            % activation functions etc.
            f = @(v) bsxfun(@rdivide,b+v, 1+sum(v,1)).^p;
            Score = @(y) -sum(sum(y)) + beta*sum(sum(y,2).^2)/2;
            Penscore = @(w) Score(f(w*X)) + alpha*sum(w(:).^2)/2; 
            Opt = @(w) OptFn(w,X,p,b,alpha, beta);
            Options = optimoptions('fmincon', 'Algorithm', 'sqp', 'GradObj', 'on', 'display', 'notify');%, 'Display', 'iter');%, 'DerivativeCheck', 'on');

            % now do multiple starts of optimization, from each gene
            Scores = zeros(nStarts,1);
            wStore = zeros(d,VectorLength,nStarts);
            
            parfor i=1:nStarts
%                 w0=zeros(d,nBims);
%                 w0(1,i) = 1;
                w0 = accumarray([1 i], 1, [d VectorLength]);
                w = fmincon(Opt, w0, [], [], [], [], zeros(size(w0)), [], [], Options);
            
                Scores(i) = Penscore(w);
                
                wStore(:,:,i) = w;
                
                if ~ScoreOnly
                    fprintf('%s score %f\n', gBim.IDsToNames(i), Scores(i));
                end
            end
            
            % this is a bit of a hack to make it return only the best score
            if ScoreOnly
                BestScore = min(Scores);
                Clu = []; Suffix1 = []; Suffix2 = []; Decision = [];
                return;
            end
                
            
            % now go through them in order, and see if user wants to split
            [ScoreSorted, ScoreOrder] = sort(Scores, 'ascend');
            LastScore = inf;
            while 1
                for i=1:VectorLength
                    seed = ScoreOrder(i);

                    % only consider this if different score to previous (to
                    % avoid seeing multiple of the same local maximum)
                    if i>1 & abs(ScoreSorted(i)-LastScore)<SkipThresh; 
                        fprintf('Skipping %s, Score %f\n', gBim.IDsToNames(seed), ScoreSorted(i));
                        continue; 
                    end
                    LastScore = ScoreSorted(i);

                    
                    w = wStore(:,:,seed);
                    
                    y = f(w*X);
                    Cells1 = find(max(y,[],1)>Thresh & y(1,:)>=y(2,:));
                    Cells2 = find(max(y,[],1)>Thresh & y(1,:)<y(2,:));
                    Cells0 = find(max(y,[],1)<=Thresh);
                    n1 = length(Cells1);
                    n2 = length(Cells2);
                    n0 = length(Cells0);
                    
                    if isempty(Cells0)
                        ColorOrder = 'br';
                    else
                        ColorOrder = 'kbr';
                    end

                    Clu = zeros(g.nCells,1);
                    Clu(Cells1) = 1;
                    Clu(Cells2) = 2;


                    % plot weights of all genes
                    figure(1); clf; subplot(2,1,1);
                    set(gcf, 'Name', PlotTitle );
                    plot(1:VectorLength, w(1,:), 'b.');
                    plot(1:VectorLength, w(2,:), 'r.');
                    set(text(1:VectorLength, w(1,:), gBim.GeneName), 'color', 'b', 'Interpreter', 'none');
                    set(text(1:VectorLength, w(2,:), gBim.GeneName), 'color', 'r', 'Interpreter', 'none');
                    set(gcf, 'Name', PlotTitle);
                    xlabel('Bimodality order');
                    ylabel('Weight');
                    ylim([0 max(w(:))]*1.1);

                    subplot(2,1,2); cla; hold on
                    wX = w*X;
                    plot(wX(1,Cells1) + rand(1,n1)*Jitter, wX(2,Cells1) +rand(1,n1)*Jitter, 'b.');
                    plot(wX(1,Cells2) + rand(1,n2)*Jitter, wX(2,Cells2) +rand(1,n2)*Jitter, 'r.');
                    plot(wX(1,Cells0) + rand(1,n0)*Jitter, wX(2,Cells0) +rand(1,n0)*Jitter, 'k.');
                    xlabel('Team 1 sum');
                    ylabel('Team 2 sum')
                    
                    figure(2); clf
                    % plotmatrix of genes with scores more than 0.1
                    [Wtsorted1 WtOrder1] = sort(w(1,:),'descend');
                    nShow1 = min(13,max(find(Wtsorted1>.1)));
                    [Wtsorted2, WtOrder2] = sort(w(2,:),'descend');
                    nShow2 = min(13,max(find(Wtsorted2>.1)));
                    if nShow1>0 & nShow2>0
                        gBim.PlotMatrix({WtOrder1(1:nShow1), WtOrder2(1:nShow2)},Clu,0,ColorOrder,5);
                        set(gcf, 'Name', [PlotTitle ': Splitting Genes']);
                    end
                    
                    % now with colors determined by BackSpin
%                     figure(3); clf
%                     gBim.PlotMatrix({WtOrder1(1:nShow1), WtOrder2(1:nShow2)},gBim.CellInfo.level2class,0,colorcube(16),'.',5);
%                     set(gcf, 'Name', PlotTitle);
%                     
%                     figure(4); clf
%                     GeneOrder = BestDiscriminants(g,Clu);
%                     g.PlotMatrix(GeneOrder(1:8),Clu,0,ColorOrder,'.',5);
%                     set(gcf, 'Name', [PlotTitle ': Best Discriminants']);

                    % show genes with best Choice probability in group 1 vs
                    % group 2
                    figure(4); clf
                    set(gcf, 'Name', PlotTitle );
                    cp21 = g.ChoiceProb(Cells2, Cells1);
                    Color = abs(cp21-.5).^pow * 2^(pow-1) .*sign(cp21-.5) + .5;

                    [~, Cporder] = sort(abs(cp21-.5), 'descend');
                    TopGenes = Cporder(1:nTops);

                    subplot(3,2,1);
                    ColoredWords(g.GeneName(TopGenes), Color(TopGenes), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 2 vs Group 1');
                    
                    Favs = g.NamesToIDs(intersect(g.GeneName,FavoriteGenes));
                    subplot(3,2,2);
                    ColoredWords(g.GeneName(Favs), Color(Favs), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 2 vs Group 1');
                    
                    % now find best choice prob of group 1 vs the rest
                    cp1 = gAll.ChoiceProb(TheseCells(Cells1));
                    Color = abs(cp1-.5).^pow * 2^(pow-1) .*sign(cp1-.5) + .5;
                    [~, Cporder1] = sort(abs(cp1-.5), 'descend');
                    TopGenes1 = Cporder1(1:nTops);
                    
                    subplot(3,2,3);
                    ColoredWords(g.GeneName(TopGenes1), Color(TopGenes1), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 1 vs the rest');
                    
                    subplot(3,2,4);
                    ColoredWords(g.GeneName(Favs), Color(Favs), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 1 vs the rest');

                    % now find best choice prob of group 2 vs the rest
                    cp2 = gAll.ChoiceProb(TheseCells(Cells2));
                    Color = abs(cp2-.5).^pow * 2^(pow-1) .*sign(cp2-.5) + .5;
                    [~, Cporder2] = sort(abs(cp2-.5), 'descend');
                    TopGenes2 = Cporder2(1:nTops);
                    
                    subplot(3,2,5);
                    ColoredWords(g.GeneName(TopGenes2), Color(TopGenes2), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 2 vs the rest');
                    
                    subplot(3,2,6);
                    ColoredWords(g.GeneName(Favs), Color(Favs), 6, BlueWhiteRed(200));
                    xlim([.4, 6.6]);
                    title('Group 2 vs the rest');
                    


                    if ~isempty(DisplayGenes)
                        figure(5); clf
                        g.PlotMatrix(DisplayGenes,Clu,0,ColorOrder,'.',5);
                        set(gcf, 'Name', [PlotTitle ': Old Favorite Genes']);
                    end
                    
%                     figure(3); clf
                    
                    
%                     xmax=max(w(1,:)*X);
%                     plot([1 xmax], [1 xmax]*(Theta^-1 - 1))
    %                 for j=1:d
    %                     Names = cellstr(gBim.IDsToNames(find(w(j,:)>.1)));
    %                     for k=1:length(Names) fprintf('%s ', Names{k}); end;
    %                     if j<d; fprintf('\n vs. '); else fprintf('\n'); end;
    %                 end
                    fprintf('Seed %s, Score %f = %f + %f + %f\n',  gBim.IDsToNames(seed), ScoreSorted(i), -sum(sum(y)), beta*sum(sum(y,2).^2)/2, alpha*sum(w(:).^2)/2); 
                    
                    Response = input('(s)plit; (a)ll again; (q)uit don''t split; return to see next: ', 's');
                    if Response=='k'
                        keyboard;
                        continue;
                    end
                    if length(Response)>0 break; end
                end
                if Response=='s' | Response=='q' | Response=='2' | Response=='1' ; break; end
            end
            
            if Response=='q' 
                Clu = [];
            end


if 0            
            figure(2); print('-dpdf', ['PDFs\Scatter ' PlotTitle '.pdf']);
            figure(4); print('-dpdf', ['PDFs\ColorGenes ' PlotTitle '.pdf']);
end
            
            % choose subgroup names
            [~,Order1] = sort(cp1, 'descend');
            [~,Order2] = sort(cp2, 'descend');
            Pos = find(Order1~=Order2, 1, 'first');
            if isempty(Pos); Pos=1; end; % if Order1 and Order2 are completely equal, give up
            NameGene1 = Order1(Pos);
            NameGene2 = Order2(Pos);

            
            if Response=='2'
                Suffix1 = g.IDsToNames(NameGene2);
                Suffix2 = g.IDsToNames(NameGene1);
                Clu=(3-Clu).*(Clu>0); % swap clusters 1 and 2, leave 0
            else                
                Suffix1 = g.IDsToNames(NameGene1);
                Suffix2 = g.IDsToNames(NameGene2);
            end
            
            Decision.w = nan(2, gAll.nGenes);
            Decision.w(:, BimOrder(1:VectorLength)) = w;
            Decision.y = nan(2, g.nCells);
            Decision.y(:, TheseCells) = y;
            Decision.Score = ScoreSorted(i);
            Decision.UserResponse = Response;
        end

        
        function Scores = SplitNull(g, n)
            % Scores = SplitNull(g, n)
            %
            % return a null distribution of best scores from splitting
            % n random permutations of a geneset g
            
            Scores = zeros(n,1);
            for i=1:n
                [~,~,~,~,Scores(i)] = Split(g.Randomize(), 1:g.nCells, [], [], 1);
            end
        end
        
        function [w, Objective, Errors] = L1Classify(g, wMax, C1, C0);
            % w = L1Classify(g, wMax, C1, C0);
            % 
            % runs a L1-constrained SVM (with no offsets and a couple of
            % other tweaks) to separate cells in class C1 from those in class C0
            % lambda: maximum 1-norm of weights (try 1)
            % 
            % C1 and C0 are targets 
            % if vectors of length nCells, then they are weights from 0 to 1
            % otherwise parsed by IdentifyCells
            % 
            % returns: w is weight vector (nGenes dimensional)
            % Errors is how much they missed the criterion by (nCells
            % dimensional)
            
            nG = g.nGenes;
            nC = g.nCells;

            if length(C1) ~= nC
                C1 = ismember(1:nC,g.IdentifyCells(C1));
            end
                
            if nargin<4
                C0 = 1-C1;
            elseif length(C0) ~= nC
                C0 = ismember(1:nC,g.IdentifyCells(C0));
            end
                        
            d = nG + nC; % number of dimensions (including slack vars)

            % Linear programming: variables x are gene weights then slack vars
            % penalties are f.x
            f = [zeros(nG,1) ; ones(nC,1)];

            % inequality constaints for each cell, plus one for weight sum
            A = zeros(nC+1,d);
            b = zeros(nC+1,1);
            
            % for cells in C1, -w.x_i - s_i < -1
            A(C1,1:nG) = -g.GeneExp(:,C1)';
            b(C1) = -1;
            % for cells in C0, w.x_i - s_i < 0
            A(C0, 1:nG) = g.GeneExp(:,C0)';
            b(C0) = 0;
            % now put in slack variables term
            A(1:nC,nG+1:nG+nC) = -eye(nC);
            % and now weight sum constraint
            A(end,1:nG) = 1;
            b(end) = wMax;

            % lower bound
            lb = zeros(d,1);

            opts = optimoptions('linprog', 'Algorithm', 'dual-simplex');
            [x, Objective] = linprog(f,A,b,[],[],lb,[],[],opts);

            w = x(1:nG);
            Errors = x(nG+1:end);
        end
            
        
        function [GeneOrder Fstat] = BestDiscriminants(g, Class, Regularize);
            % [GeneOrder Fstat] = BestDiscriminants(g, Class);
            % sorts the GeneIDs in order of how well they separate the
            % cells into classes. (default is g.Class)
            % 
            % uses ANOVA F-statistic (between group var/within group var)
            %
            % optional input Regularize is added to the denominator of the
            % F statistic
             
            if nargin<2
                Class = g.Class;
            end
            
            if nargin<3
                Regularize = 1;
            end
            
            [c, ia, ic] = unique(Class);
            nGroups = length(c);
            
            Means = zeros(g.nGenes,nGroups); 
            Resids = zeros(g.nGenes, g.nCells);
            nCells = zeros(nGroups,1);
            for i=1:nGroups
                % find mean of each group, and residuals of all matrix elements
                GpCells = find(ic==i);
                Means(:,i) = mean(g.GeneExp(:,GpCells),2);
                Resids(:,GpCells) = bsxfun(@minus,g.GeneExp(:,GpCells),Means(:,i));
                nCells(i) = length(GpCells); % number of cells in each group
            end
            
            UnexplainedVar = sum(Resids.^2,2)/(g.nCells-nGroups);
            OverallMean = mean(g.GeneExp,2);
            ExplainedVar = bsxfun(@minus,Means,OverallMean).^2 * nCells / (nGroups-1);
            
            Fstat = ExplainedVar./(UnexplainedVar+Regularize);
            
            MaxMean = max(Means,[],2);
            MarkovStat = MaxMean./mean(g.GeneExp,2);
            [~, GeneOrder] = sort(Fstat,'descend');
%            [sorted, GeneOrder] = sort(MarkovStat,'descend');
            
            % now try largest group mean over total mean (only works for two groups)
%             for i=1:20
%                 gg = GeneOrder(i);
%                 fprintf('%s: F=%.3f, Markov=%.3f\n', g.IDsToNames(gg), Fstat(gg), MarkovStat(gg));
%             end
        end
        
        function [cp Order] = ChoiceProb(g, C1,C2)
            % [cp Order] = ChoiceProb(C1,C2)
            %
            % Computes each gene's choice probability for cell class C1 vs
            % C2. i.e. the probability that a random cell in class C1 will
            % express more of each gene than a random cell in class C2.
            % (Uses the rank sum statistic)
            % 
            % C1 and C2 can be logical, or a set of indices, or a string in
            % which case it specifies an entry in Class to be
            % compared using strmatch.
            %
            % Default for C2 is everyone not in C1. If C1 and C2 overlap,
            % duplicates are removed from C2 (not from C1).d
            % 
            % return value cp is between 0 and 1 for all genes
            % optional output Order returns GeneIDs arranged from the most
            % C1 first. 
            
            if isstr(C1)
                C1 = strmatch(C1,g.Class);
            end
            
            if islogical(C1)
                C1 = find(C1);
            end
            
            if nargin<3
                C2 = setdiff(1:g.nCells,C1);
            end
            
            if isstr(C2)
                C2 = strmatch(C2,g.Class);
            end

            if islogical(C2)
                C2 = find(C2);
            end
            
            % remove duplicates
            C2 = setdiff(C2, C1);
            
            % concatenate
            C = [C1(:) ; C2(:)];
            
            % how many cells in each class
            n1 = length(C1);
            n2 = length(C2);
            
            
            cp = zeros(g.nGenes,1);
            parfor i=1:g.nGenes
                Ranks = tiedrank(g.GeneExp(i,C)); % because tiedrank works on columns
                Ranksum = sum(Ranks(1:n1));
            
                cp(i) = (Ranksum - n1*(1+n1)/2)/n1/n2; % formula for choice prob
            end
            
            if nargout>=2
                [~,Order] = sort(cp, 'descend');
            end
                       
        end
        
        function gSub = SubHierarchy(g, BaseClass)
            % gSub = SubHierarchy(BaseName)
            %
            % Finds all cells whose Class is a subclass of BaseClass 
            % (i.e. starts with it). e.g. BaseClass = '.Cnr1.Reln', no
            % trailing period.
            % 
            % Then strips BaseClass out of their 
            % ClassName - unlike CellSubset, which doesn't.
            
            MyCells = strmatch(BaseClass, g.Class);
            
            gSub = g.CellSubset(MyCells);
            
            for i =1:gSub.nCells
                gSub.Class{i} = gSub.Class{i}(length(BaseClass)+1:end);
            end
            
        end
        
        function gTrunc = TruncateHierarchy(g, nLevels)
            % gTrunc = TruncateHierarchy(nLevels)
            %
            % Truncates a hierarchy to at most nLevels levels deep
            
            gTrunc = g;
            for i=1:g.nCells
                MyClass = g.Class{i};
                Points = find(MyClass=='.', nLevels+1);
                if length(Points)<nLevels+1 % not that many levels for this cell
                    gTrunc.Class{i} = MyClass;
                else
                    gTrunc.Class{i} = MyClass(1:max(Points)-1);
                end
            end
        end
                    
                
        
        function TopGenes = HierarchyView(g, StartClass, LastLevel, pow, GenesToShow, GenesToKeep);
            % TopGenes = HierarchyView(StartClass, LastLevel, pow, Genes2Show)
            % 
            % for a hierarchical classification in g.Class
            % makes a dendrogram of gene names in red and blue for showing
            % how much each subclass expresses each gene relative to the
            % whole population. 
            %
            % Starts analyzing the hierarchy from StartClass (default '' -
            % make sure you don't include a trailing '.')
            %
            % Goes up to LastLevel level of hierarchy (default all). If
            % this is a length-2 vector, the second element says where to
            % stop shrinking the vertical gap between boxes (otherwise they are equal)
            %
            % pow is power to raise the choice prob statistic to when
            % chosing a color
            %
            % GenesToShow is which genes to display - note you
            % can define a FavoriteGenes.m, in which case that is default.
            %
            % GenesToKeep is how many top per division to collect for the
            % output, TopGenes. For this to work, GenesToShow must be {}
            %
            
            % visual parameters
            Landscape = 0;
            FontSize = 9;
            nColumns = 13; % number of columns of text within each box
            
            
            if nargin<2;
                StartClass = '';
            end
            FirstLevel = sum(StartClass=='.');
                        
            if nargin<3 || isempty(LastLevel) || LastLevel(1)==0 || ~isfinite(LastLevel(1))
                LastLevel = max(cellfun(@(x) sum(x=='.'), g.Class));
            end
            if length(LastLevel)==2
                StopShrinkingAt = LastLevel(2);
                LastLevel = LastLevel(1);
            else
                LastLevel = LastLevel(1);
                StopShrinkingAt = LastLevel(1);
            end
            
            if nargin<4 || isempty(pow)
                pow = 1.5;
            end
            
            if nargin<5
                GenesToShow = FavoriteGenes;
            end
            
            if nargin<6
                GenesToKeep = 10;
            end
            
            EmptyStart = (g.CellSubset(StartClass).nCells == g.nCells); % no need to show first box
            
            nLevels = LastLevel - FirstLevel + 1;
            nShrinks = min(nLevels, StopShrinkingAt+1);
            
            Bottom = -0/8; % edges of the last boxes in the right column
            Top = 1.01;
            xGap = .3; % how much to shift all boxes right, except last row which is shifted 1

            % box size scaling parameters
            TotLength = 2  + (nLevels -2 -EmptyStart)*xGap;
            Height = 1/2^(nShrinks-1)/(Top-Bottom) - .025; % of single box
            Width = 0.95/TotLength; % of single box
            xSpace = xGap/TotLength; % distance between boxes except last col
            LastxSpace = 1/TotLength; % distance between boxes
            xShift = (0.025-xGap*EmptyStart)/TotLength; % position of left edge of left plot
            yShift = .5/(Top-Bottom) - Bottom; %center of left plot
            ySpace = .5/(Top-Bottom); % spacing between center of second two plots
            SmallerBy = .5; % How much to shrink vertical gap with each level
            
            
            if EmptyStart
                [BinaryNames, ClassRoots] = BinarizeClassNames(g.Class);
                StartBinary = '';
            else
                [BinaryNamesPlus, RootsPlus] = BinarizeClassNames([StartClass; g.Class]);
                BinaryNames = BinaryNamesPlus(2:end);
                StartBinary = BinaryNamesPlus{1};
                ClassRoots = RootsPlus(2:end);
            end
            
            
            clf; 
            if Landscape
                set(gcf, 'papertype', 'a4', 'paperorientation', 'landscape', 'inverthardcopy', 'off', 'color', [1 1 1]);
                set(gcf,'paperunits', 'centimeters', 'PaperPosition',[.7 .7 28.5 19.5 ]);	%this is for A4
            else
                set(gcf, 'papertype', 'a4', 'paperorientation', 'portrait', 'inverthardcopy', 'off', 'color', [1 1 1]);
                set(gcf,'paperunits', 'centimeters', 'PaperPosition',[.7 .7 19.5 28.5]);	%this is for A4
            end

            
            % now loop over levels
            TopGenesAll = {};
            for Level = 0:(LastLevel-FirstLevel)
                % loop over potential branches within that level
                for i=1:2^Level;
                    MyBinary = dec2bin(i-1,Level);
                    BinaryRoot = [StartBinary, MyBinary];
                    MyCells = strmatch(BinaryRoot, BinaryNames);
%                     MyClassName = ClassRoots{Level+1}{i};
                    
                    if isempty(MyCells), continue; end;
                    
                    % now find left axis of box in the right place
                    IsRightColumn = (Level==LastLevel-FirstLevel);
                    Left = xShift + xSpace*Level + (LastxSpace-xSpace)*IsRightColumn;
                                    
                    if length(MyCells) == g.nCells % skip empty box
                        
%                         yCenter = yShift;
%                         axes('position', [Left, yCenter-Height/2, Width, Height]);
% %                         box on
%                         set(gca, 'xtick', []);
%                         set(gca, 'ytick', []);
%                         set(gca, 'xticklabel', '')
%                         set(gca, 'yticklabel', '')
                    else
                        PlusMinus = (MyBinary-'0')*2 - 1;
                        yCenter = yShift + ySpace*(SmallerBy.^min(1:Level, StopShrinkingAt)) * PlusMinus';
                        axes('position', [Left, yCenter-Height/2, Width, Height]);
                        
                        [h gn] = g.ColorGeneNames(MyCells, GenesToShow, pow, nColumns);
%                         title(ClassRoots{Level+FirstLevel}{i});
                        %if length(gn)>20; gn=gn(1:20); end;
                        TopGenesAll = vertcat(TopGenesAll, gn(1:GenesToKeep));
%                         TopGeneLabelsAll = vertcat(TopGeneLabelsAll, gl);
                        set(h, 'FontSize', FontSize);
                        
                    
                        % draw a rounded rectange
                        % box on
                        axis off;
                        xl = xlim; yl = ylim;
                        rectangle('Position', [xl(1) yl(1) diff(xl) diff(yl)], 'Curvature', [.05 .05], 'LineWidth', .2);

                        drawnow;
                    end
                        
                end
            end
            
            [TopGenes ia] = unique(TopGenesAll);
%             TopGeneLabels = TopGeneLabelsAll(ia);
            
        end


        function [h, GeneNames] = ColorGeneNames(g, Group, ShowGenes, pow, nColumns)
            % h = ColorGeneNames(Cells, Genes, pow, nColumns)
            %
            % shows which genes express in Cells above and below the rest 
            % of the population. (Red means more, blue less.)
            %
            % Cells is a subset of cells, can be only one cell if you like
            % if you pass a string, it uses this subhierarchy
            %
            % Genes is which genes to consider - default is to pick 100 most
            % different.
            %
            % pow is a "gamma correction" - try 2 or 3 (larger means only
            % strongest are shown)
            %
            % returns a vector of handles to each word and a list of genes
            % (in case you had asked it to find the best)

            
            if nargin>=3 & length(ShowGenes)==1
                nDefault = min([40, g.nGenes, ShowGenes]); % how many top ones to show if not specified
            else
                nDefault = min([40, g.nGenes]); % how many top ones to show if not specified
            end
            Cells = g.IdentifyCells(Group); 
            
            if nargin<3 | isempty(ShowGenes) | length(ShowGenes)==1
                RelExp = g.RelativeExpression(Cells);
                [~, order] = sort(abs(RelExp), 'descend');
                Genes = order(1:nDefault);
            else
                Genes = intersect(ShowGenes, g.GeneName);
                if length(Genes)<length(ShowGenes)
                    fprintf('ColorGeneNames: showing %d of %d genes requested\n', ...
                        length(Genes), length(ShowGenes));
                end
            end
            
            
            
            if nargin<4
                pow = 1.5;
            end
            
            if nargin<5
                nColumns = 7;
            end
            RelExp = g.RelativeExpression(Cells);
            
            RelExpMine = RelExp(g.NamesToIDs(Genes));
            GeneNames = g.IDsToNames(Genes);
            [~, PlotOrder] = sort(abs(RelExpMine), 'ascend');
            
            % compute color as number between 0 and 1
            Color = abs(RelExpMine).^pow /2 .*sign(RelExpMine) + .5;

            % plot them in the right order (white on the bottom)
            cla;
            hh = ColoredWords(GeneNames, Color, nColumns, BlueWhiteRed(200), PlotOrder);
            if nargout>=1; h = hh; end;
            xlim([.4, nColumns+.6]);
            
            
            for i=1:length(hh)
                set(hh(i), 'ButtonDownFcn', {@MDSCallback, GeneNames{i}});
            end
            
        end
        
        function ViewAllCells(g, Cells)
            % ViewAllCells(Cells)
            % iterate through all cells in the database, showing a colored
            % words plot for them. If you provide argument Cells it does
            % only those - if it is a string, it does that class.
            %
            
            if nargin<2; Cells=1:g.nCells; end
            if isstr(Cells); Cells = strmatch(Cells, g.Class); end;
            
            for c=Cells(:)'
                clf;
                
%                 subplot(1,2,1);
%                 g.ColorGeneNames(c);
%                 title('top genes');
%                 
%                 subplot(1,2,2);
                g.ColorGeneNames(c, FavoriteGenes);
                title('Favorite genes');
                
                set(gcf, 'name', sprintf('Cell %d, class %s', c, g.Class{c}));
                pause;
            end
        end   
        
        function ConfusionMatrix(g1, g2, Id1, Id2);
            % ConfusionMatrix(g2, Id1, Id2)
            % 
            % plots a confusion matrix for how cells are classified in a
            % second geneset, g2. g2 has to contain exactly the same cells, 
            % but not necessarily in the same order; Id1 and Id2 are sets
            % of unique identifiers for each cell, defaulting to g.CellName
            % and g2.CellName
            
            if nargin<3
                Id1 = g1.CellName;
            end
            if nargin<4
                Id2 = g2.CellName;
            end
            
            % sort the Ids to get in order;
            [Sorted1 ,Order1] = sort(Id1);
            [Sorted2 ,Order2] = sort(Id2);
            if ~isequal(Sorted1, Sorted2)
                error('IDs don''t match!');
            end
            
            NewOrder(Order1) = Order2; % inverse of one permutation times another
            [tab0,~,~,labels] = crosstab(g1.Class, g2.Class(NewOrder));
            tab = bsxfun(@rdivide, tab0, sum(tab0,1));
            [n1, n2] = size(tab);
            
            % decide order the second class by similarity to the first
            MeanPosition = (1:n1)*tab;% ./ sum(tab,1);
            [~, xOrder] = sort(MeanPosition);
            
            hold off
            imagesc(tab(:,xOrder)); 
            set(gca, 'xtick', 1:n2, 'xticklabel', labels(xOrder,2), 'xticklabelrotation', 90);
            set(gca, 'ytick', 1:n1, 'yticklabel', labels(1:n1,1));
            hold on
            for i=1.5:n1
                plot(xlim, [i i], 'w', 'linewidth', .05);
            end
            for i=1.5:n2
                plot([i i], ylim, 'w', 'linewidth', .05);
            end
            
            
%             grid on;
%            set(gca, 'gridcolor', [1 1 1]);
        end
        
        function b = Bimodality(g)
            % b = Bimodality();
            % finds for every gene a measure of bimodality, given by the
            % probability of seeing 0 expression, minus the value predicted
            % by a negative binomial fit from values 1 and above.
            % (see NbinMixFit)
            %
            % also insists on at least 5 cells expressing more than 4
            
            
            b = zeros(g.nGenes,1);
            parfor i=1:g.nGenes
                %if mod(i,1000)==0; fprintf('Bimodality: gene %d\n', i); end;
                MyExp = g.GeneExp(i,:);
                if sum(MyExp>4)<5 
                    b(i) = 0;
                else

                    %[r, p, s] = NbinMixFit(g.GeneExp(i,:));
%                     if isfinite(s)
%                         b(i) = 1-s;
%                     else
%                         b(i) = 0;
%                     end
                    [~,~,b(i)] = NbinFitMinP(MyExp,.2,.99);
                end
            end
            %fprintf('\n');
        end
        
        function h = ClearNames(g)
            % h = ClearNames(g)
            % gets rid of the CellName field - because IdentifyCells matches both
            % classes and names, and sometimes you only want to match
            % classes
            
            h = g;
            h.CellName = [];
        end
            
        function c = IdentifyCells(g, s, Exclude)
            % c = IdentifyCells(s, Exclude)
            % does some indexing etc to let you quickly find sets of
            % cells
            %
            % if s is a set of integers, these are the IDs
            % if s is a set of logical values, it finds these
            % if s is a string or cell array of strings, it finds any cells 
            % matching the cell name or class
            % if s is a numerical vector, it goes right through (much faster)
            %
            % the optional last argument Exclude says don't use those cells
            %
            % NB this uses union, so each cell is returned once, and they
            % might not be in the same order.
            
            % turn singleton string into 1-element cell array
            %if isstr(s); s = {s}; end
            
            if iscell(s) % cell array means union of all things in it
                c = [];
                for i=1:length(s)
                    c = [c; g.IdentifyCells(s{i})];
%                     c = union(c, strmatch(s{i}, g.Class));
%                     c = union(c, strmatch(s{i}, g.CellName));
%                      c = [c; strmatch(s{i}, g.Class)];
%                      c = [c; find(strcmp(s{i}, g.CellName))];

                end
                c = unique(c, 'stable'); % preserve order
            elseif isstr(s)
                c = [strmatch(s, g.Class) ; find(strcmp(s, g.CellName))];
            elseif islogical(s) % logical array means keep all with 1
                c = find(s); 
            elseif isnumeric(s) % numeric array means those cells
                c = s;
            else
                error('IdentifyCells: not sure what to do with input');
            end
            
            % remove any excluded cells
            if nargin>=3
                c2 = g.IdentifyCells(Exclude);
                c = setdiff(c, c2);
            end
            
        end
            

                
        function [d, m, p] = RelativeExpression(g, Group1, Group2, reg)
            % [d, m, p] = RelativeExpression(g, Group1, Group2, reg)
            %
            % computes a relative expression index for cells of Group1 vs
            % cells of Group2. 
            % 
            % Group1 and Group2 inputs are parsed by IdentifyCells() - but
            % if they are numeric vectors of length nGenes, those are used
            % as data directly.
            %
            % Default for Group2 is all cells not in Group1.
            %
            % Output d is (Mean1 - Mean2)./(Mean1 + Mean2 + reg) for every
            % gene. Default value for reg is 1. Mean1 and Mean2 are trimmed
            % means (34%) of each cell in each group. NOW MEDIANS
            % TEMPORARILY
            %
            % Output m is (Mean1+Mean2)/2
            %
            % optional third output is the p-value of a rank sum test
            % asking whether each gene is expressed differently between
            % groups
       
%             if length(Group1)==g.nCells
%                 Mean1 = Group1(:);
%             else
                Group1=g.IdentifyCells(Group1);
                Mean1 = median(g.GeneExp(:,Group1)');
                %Mean1 = trimmean(g.GeneExp(:,Group1)',34, 'round', 1);% + reg;
%             end
            
            if nargin<3 || isempty(Group2)
                Group2 = setdiff(1:g.nCells, Group1);
                Mean2 = median(g.GeneExp(:,Group2)');
                %Mean2 = trimmean(g.GeneExp(:,Group2)',34, 'round', 1);% + reg;
            elseif length(Group2)==g.nCells
                Mean2  = Group2(:);
            else
                Group2=g.IdentifyCells(Group2);
                Mean2 = median(g.GeneExp(:,Group2)');
                Mean2 = trimmean(g.GeneExp(:,Group2)',34, 'round', 1);% + reg;
            end

            if nargin<4 || isempty(reg); reg=1; end
            
            
%            m = (Mean1 + Mean2)/2;
            m = max(Mean1,Mean2);
%            m = trimmean(g.GeneExp',34);
            d = (Mean1-Mean2)./(Mean1+Mean2+reg);
            
            if nargout>=3
                p = zeros(g.nGenes,1);
                for i=1:g.nGenes
                    p(i) = ranksum(g.GeneExp(i,Group1), g.GeneExp(i,Group2));
                end
            end
        end
            
        function [h PlotMe] = DichotomyPlot(g, Group1, Group2, reg);
            % [h, Plotted] = DichotomyPlot(Group1, Group2, reg);
            % analyzes genes that differ particular between Groups (see
            % RelativeExpression for more info on input format)
            % 
            % only plots those differing at p<.05
            %
            % returns handles of names plotted, plus indices of all plotted genes
                        
            if nargin<3; Group2=[]; end
            if nargin<4; reg=0; end
            [d, m, p] = g.RelativeExpression(Group1, Group2, reg);
            
            PlotMe = find(p<.05);
            %PlotMe = 1:g.nGenes;            
            
            semilogx(m(PlotMe),d(PlotMe), 'k.');
            ylim([-1 1]);

            h = text(m(PlotMe)*1.05, d(PlotMe), g.GeneName(PlotMe), 'Interpreter', 'none');
if 1            
            for i=1:length(h)
                set(h(i), 'ButtonDownFcn', {@MDSCallback, g.GeneName{PlotMe(i)}});
            end
end
            set(h, 'color', 'b');
            grid on
            
            xlabel('Mean Expression');
            ylabel('Relative Expression Index');
            grid on;

            if nargout==0; clear h; end

        end
        
        function h = Binarize(g, Threshold)
            % returns a binarized version, with 1 if expression>=Threshold,
            % 0 if expression<=Threshold
            %
            % Default Threshold is 5.
           
            if nargin<2; Threshold = 5; end
            h = g;
            h.GeneExp = double(g.GeneExp>=Threshold);
        end
        
        function h = ComputetSNE(g, WhichGenes, nDim);
            % h = tSNE(nDim);
            % does a nDim-dimensional tSNE plot using the external drtoolbox
            % and saves the output in the tSNE variable
            % default dim = 2
            
            if nargin<3
                nDim =2;
            end
            
            Reg = 0.1;
            r = 2;
            
            x = g.Exp(WhichGenes,1:g.nCells) + Reg; % expression of active genes: nActive by nC
            p = x./(x + r); % negbin parameter:  nActive by nC
            L = bsxfun(@plus, x'*log(p), sum(r*log(1-p), 1)); % L(i,j): log likelihood of cell i in cluster defined by cell j

%            LogPji = bsxfun(@minus, L, LogSumExp(L)); % log conditional probability of cluster j for cell i
             LogPij = bsxfun(@minus, L', LogSumExp(L')); % log conditional probability of cluster i for cell j
%             LogPsym = LogPji + LogPij;
            LogPsym = LogPij + LogPij';
            
            [P, beta] = d2p(LogPsym/median(LogPsym(:)));
            
            % make Rodriguez-Liao plot
            Dist = ClosestBigger(-LogPsym, beta);
            Size = 1*ones(size(Dist));
            [gp, gpn] = grp2idx(g.Class);
            for i=1:length(gpn)
                Mine = find(gp==i);
                [~, TopIdx] = max(beta(Mine));
                Size(Mine(TopIdx)) = 25;
                % find closest point in another cluster
            end
            
            figure(2389074); clf;
            g.BrushableScatter([(beta), (Dist)], Size);
            grid on
            xlabel('density');
            ylabel('Distance to closest denser point');
            title('Rodriguez-Laio plot');
            
            
            
            figure(2389075); clf;
            [~, ~, ClassNum] = unique(g.Class, 'stable');
            InitialSolution = [zscore(ClassNum(:)), randn(g.nCells,nDim-1)*1e-4];
            th = 2*pi*(ClassNum(:) - min(ClassNum))/(max(ClassNum) - min(ClassNum));
            InitialSolution = [cos(th), sin(th)];
            h = g;
            h.tSNE = tsne_p(P, ClassNum, InitialSolution)';
            
            

%             h = g;
%             h.tSNE = compute_mapping(g.GeneExp', 'tSNE', nDim)';
        end
        
        function BrushableScatter(g, x, Size, Sym, Clr, Size0)
            % h = BrushableScatter(g, x, Size, Sym, Clr, Size0)
            %
            % makes a brushable 2 or 3d scatterplot from data matrix x (nCells by nDims). 
            % Use WhosBrushed to find out the names of brushed cells
            % colors and glyps chosen according to class. 
            % optional input Size (default 3) can be scalar or vector of
            % size nPoints - if 0, always plot points
            %
            % If Sym and Clr are provided, then they are arrays to look up
            % the appropriate color and size for each symbol - use
            % containers.map
            % 
            % Size0 says what to do with points of Size==0; if Size0 is 0,
            % they won't be plotted, otherwise they are points of this size
            
            
            if nargin<3; Size = 10; end
            if nargin<6; Size0= 2; end
            
            if length(Size)==1
                Size = Size * ones(g.nCells, 1);
            end
            
            if iscell(x) & numel(x)==2
                % passing genenames
                x = [g.Exp(x{1})', g.Exp(x{2})'];
            end
            
            Group = g.Class;
            if isempty(Group) 
                Group = ones(g.nCells,1);
            end
            
            [uGrpNm, ~, uGrpInd] = unique(Group, 'stable');
            nGroups = length(uGrpNm);
            
            if nargin<4 | isempty(Sym) | isempty(Clr)
                % make symbols alternating colors and shapes
                Sym0 = cellstr(repmat('o+*hxsd^<v>p', [1 ceil(nGroups/12)])');
                Clr0 = num2cell(HsvNotYellow(ceil(nGroups*1.2)), 2); % turn down green so can see yellow, 1.2 so ends not identical
                
                Sym = containers.Map(uGrpNm, Sym0(1:nGroups));
                Clr = containers.Map(uGrpNm, Clr0(1:nGroups));
            end

            % plot it one group at a time
            holding = ishold;
            if ~holding
                 cla; 
                 hold on; 
            end
            set(gca, 'clipping', 'off'); axis off
%             cla; hold on
            ThreeD = (size(x,2)>=3);
            
%             h = zeros(nGroups,1);
            for i=1:nGroups
                gn = uGrpNm{i};
                pts = find(uGrpInd(:)==i & Size(:)>0);
                pts0 = find(uGrpInd(:)==i & Size(:)<=0);
                MySize = Size(pts);
                
                if ThreeD
%                        h = plot3(x(pts,1), x(pts, 2), x(pts, 3), Sym(i), 'Color', Clr(i,:));
                        h(i) = plot3(x(pts,1), x(pts, 2), x(pts, 3), Sym(gn), 'Color', Clr(gn));
                else
                        h(i) = scatter(x(pts,1), x(pts,2), MySize, Clr(gn), Sym(gn));
                        if Size0>0
                            hh(i) = scatter(x(pts0,1), x(pts0,2), Size0, Clr(gn), '.');
                            set(hh(i), 'UserData', g.CellName(pts0));
                        end
%                     else
%                         h = plot(g.tSNE(1,pts), g.tSNE(2,pts), Sym(i), 'Color', Clr(i,:));
                end
                set(h(i), 'UserData', g.CellName(pts));
            end
            
            if nGroups>1
                legend(flipud(h(:)), flipud(uGrpNm(:)));
            end
            
            if ~holding
                hold off
            end
        end
        
        function tSNEPlot(g, Gene, SizeMult, SizeOffset, LogAdd, Sym, Clr, Size0)
            % tSNEPlot(Gene, SizeMult, SizeOffset, LogAdd, Sym, Clr, Size0)
            % plots the tSNE analysis - first call ComputetSNE
            % 
            % optional argument plots marker size as 
            % SizeOffset+SizeMult*(log(LogAdd+ Exp) - log(LogAdd))
            % zeros plotted as Size0 (default 0 means not plotted)
            % if Gene is cell array: will plot all of them in turn
            % if it is blank, plot all at size Size0
            %
            % Sym and Clr optional arguments passed to BrushableScatter
            
            if nargin<2; Gene = []; end
            if nargin<3 | isempty(SizeMult); SizeMult = 15; end
            if nargin<4 | isempty(SizeOffset); SizeOffset = 0; end
            if nargin<5 | isempty(LogAdd); LogAdd = 2; end
            if nargin<6; Sym = []; end
            if nargin<7; Clr = []; end;
            if nargin<8 | isempty(Size0); Size0 = SizeMult; end;
            
            if isstr(Gene) & ~isempty(Gene) & Gene(end)=='*'
                Gene = g.NameStartsWith(Gene(1:end-1));
            end
            
            if iscell(Gene) 
                for i=1:length(Gene)
                    cla;
                    g.tSNEPlot(Gene{i},SizeMult, SizeOffset, LogAdd, Sym, Clr, Size0);
                    pause;
                end
                return
            end
            
            if isnumeric(Gene) && length(Gene)==g.nCells % array goes straight thru
                Exp = Gene;
                Size = SizeOffset+SizeMult*(log(1+Exp/LogAdd));
            elseif isnumeric(Gene) & ~isempty(Gene) % convert IDs to Names
                g.tSNEPlot(g.IDsToNames(Gene),SizeMult, SizeOffset, LogAdd, Sym, Clr, Size0);
                return
            elseif ischar(Gene) & ~isempty(Gene)
                Exp = g.Exp(Gene);
                Size = SizeOffset+SizeMult*(log(1+Exp/LogAdd));
                if isempty(Size)
                    warning(sprintf('Gene %s not found', Gene));
                    return
                end
                
            elseif isempty(Gene)
                if nargin<8 
                    Size = 20;
                else
                    Size = Size0;
                end
            else
                error('Not sure what to do with input!\n');
            end
            
            %cla; 
            hold on; set(gca, 'clipping', 'off'); axis off
            g.BrushableScatter(g.tSNE', Size, Sym, Clr, Size0);

            if nargin>=2 && (ischar(Gene) || isscalar(Gene)) && ~isempty(Gene)
                title(sprintf('Size: %s expression', g.IDsToNames(Gene)));
%             else
%                 title('tSNE plot');
            end
        end
        
        
        function g = MDSPlot(g, CenterGene, nShowMax, p, ShowAnticorrel, nDims, MinCorrShow)
            % g = MDSPlot(CenterGene, nShowMax, p, ShowAnticorrel, nDims)
            %
            % do a multidimensional scaling plot showing gene relationships
            % among the nShow (default 30) most correlated with CenterGene
            %
            % Click on a gene and it takes you to OMIM and Allen !
            %
            % similarity measure is max(Correlation,0).^p (default p=4)
            % 
            % nDims can be 2 or 3 - default 2
            %
            % lines join top 10% of correlations amongst this set
            %
            % output argument is because it computes a correlation matrix -
            % if you call this multiple times do g = g.MDSPlot(...)
            
            if nargin<3, nShowMax =30; end
            if nargin<4, p=4; end
            if nargin<5, ShowAnticorrel = 1; end
            if nargin<6, nDims = 2; end
            if nargin<7, MinCorrShow = 0; end;
            
            CenterID = g.NamesToIDs(CenterGene);
            
            % compute correlation between each gene
            if isempty(g.GeneCorrels)
%                 gEx = g.GeneExp';
%                 g.GeneCorrels = corrcoef(gEx);
                z = zscore(g.GeneExp,1,2);
                g.GeneCorrels = z*z'/g.nCells;
                g.GeneCorrels(~isfinite(g.GeneCorrels))=0;
            end
            
            % compute similarity measure
            
             sim = min(max(g.GeneCorrels,0).^p,1);
%             nsim = sqrt(sum(sim.^2));
%             sim2 = sim^2./(nsim'*nsim);
%             sim = sim2; sim(1:(1+g.nGenes):g.nGenes^2)=1;
            
            
            % find nShow best related genes
            [sorted, order] = sort(g.GeneCorrels(:,CenterID), 'descend');
            nShow = min(nShowMax, sum(sorted>=MinCorrShow));
    
            if ShowAnticorrel; clf; end
            for p=1:(1+ShowAnticorrel)
                if ShowAnticorrel; subplot(1,2,p); end
                if p==1
                    Friends = order(1:nShow);
                else
                	Friends = order(g.nGenes-nShow:g.nGenes);
                    set(gca, 'Color', [0.3 0.3 0.3]);
                end
                

                % subset of similarity matrix
                MySim = sqrt(1-sim(Friends, Friends)); % change to dissimilarity
                MySim(1:1+size(MySim,1):end)=0; % make diags exactly 0 or it will crash
                MySim = (MySim+MySim')/2; % make symmetric or it will crash

                % do the mds calculation
                try
                    Y = mdscale(MySim,nDims);
                
                    % plot points for each gene in 2 or 3 dims
                    hold on
                    if nDims==2
                        plot(Y(:,1), Y(:,2), 'k.');    
                    else
                        plot3(Y(:,1), Y(:,2),  Y(:,3),'k.');
                    end
                    axis tight;

                    % colormap for showing top 10% of correlations
                    cmap = flipud(hot(100));

                    % find 90th and 100th percentile similarity
                    MyCorrels = g.GeneCorrels(Friends,Friends); 
    %                MyCorrels(1:(nShowMax+1):nShowMax^2)=0;
                    MyCorrels(1:(nShow+1):nShow^2)=0;
                    pCor = prctile(MyCorrels(:), [90 100]); 
                    ShowMinCorr = pCor(1);
                    ShowMaxCorr = pCor(2);
                    % ShowMinCorr = .7;
                    % ShowMaxCorr = .8;
    %                 ShowMinCorr = MinCorrShow;
    %                 ShowMaxCorr = .8;

                    % now plot lines
                    for Matches=1:length(Friends)
                        for k=1:length(Friends)
                            if Matches==k; continue; end;

                            ThisCorrel = MyCorrels(Matches, k);

                            if ThisCorrel>ShowMinCorr
                                cIndex = ceil(100*min((ThisCorrel-ShowMinCorr)/(ShowMaxCorr-ShowMinCorr),1));
                                Col = cmap(cIndex,:);
                                if nDims==2
                                    plot(Y([Matches k], 1), Y([Matches k], 2), 'color', Col);
                                else
                                    plot3(Y([Matches k], 1), Y([Matches k], 2), Y([Matches k], 3), 'color', Col);
                                end
                            end
                        end
                    end

                    % now write gene names on 
                    if nDims==2
        %                h = text(Y(2:end,1), Y(2:end,2), g.GeneName(Friends(2:end)));
                        h = text(Y(:,1), Y(:,2), g.GeneName(Friends), 'Interpreter', 'none');
        %                 text(Y(1,1), Y(1,2), g.GeneName(CenterID), 'color', 'b');
                    else
                        h = text(Y(:,1), Y(:,2), Y(:,3), g.GeneName(Friends), 'Interpreter', 'none');
        %                 text(Y(1,1), Y(1,2), Y(1,3), CenterGene, 'color', 'b');
                    end

                    % color by correlation - center gene in blue
                    nColors = 100;
                    cMap = flipud(jet(nColors+1));
                    for i=1:length(Friends)
                        color = cMap(floor((g.GeneCorrels(Friends(i),CenterID)+1)/2*nColors) + 1,:);
                        if ShowAnticorrel; set(h(i), 'color', color); end
                    end
    %                 if p==1;
    %                     set(h(1), 'color', 'b');
    %                 end

                    for i=1:length(h)
                        set(h(i), 'ButtonDownFcn', {@MDSCallback, g.GeneName{Friends(i)}});
                    end

                    title(sprintf('Lines show correlations of %f to %f', ShowMinCorr, ShowMaxCorr));
        %            set(gca, 'ButtonDownFcn', {@MDSCallback, Y, g.GeneName});

                    set(gca, 'xtick', []);
                    set(gca, 'ytick', []);
                    box on
                catch me
                    fprintf('Error doing mds: %s\n', me.getReport);
                end

            end
        end

    end
    
end

function MDSCallback(src, evt, text)
    system(['start chrome "http://www.ncbi.nlm.nih.gov/omim/?term=' text '"']);
    system(['start chrome "http://mouse.brain-map.org/search/show?search_type=gene&search_term=' text '"']);
end
       
