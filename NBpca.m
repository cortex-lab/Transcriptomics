function [w, x, MargLike] = NBpca(y, r, nF, x0, xStart, xPrior, wPrior);
% [w, x, MargLike] = NBpca(y, r, nF, x0, xStart, xPrior, wPrior);
%
% fits a negative binomial PCA model to a nObs by nDims dataset y
%
% p = 1./(1+exp(x*w)), so predicted mean is r*exp(w*x) 
%
% output x is a nObs by nF matrix. nF is number of adjustable factors
% input x0 is fixed factors, nObs by nF0. By default a single column of all ones
% but it could be anything (such as class identifiers). 
%
% output w is a (nF+nF0) by nDims matrix of weights, with the last rows being
% biases or other things that multiply the fixed factors
%
% output MargLike is marginal likelihood of all points under the model (for the
% best values of w, not integrated over w). MAYBE INCORRECT, CHECK THIS
%
% optional input xStart is start point (nObs by nF), default is to use PCA
%
% optional inputs xPrior and wPrior are the value for the diagonal of a 
% concentration matrix regularize the two fits (but not fixed/bias terms)
% - default is 0 for x, in which case the x variables are instead svd'd 
% and 1e-4 for w, just numerical regularization.
%

MaxIter = 100;
Tol = .1;
Verbose = 1;

[nObs, nDims] = size(y);

if nargin<4 || isempty(x0)
    x0 = ones(nObs,1);
    nF0 = 1;
else 
    nF0 = size(x0, 2);
end

if nargin<5 || isempty(xStart)
% initialize by PCA of log data if zero
    z0 = log(y + r);
    z1 = z0 - x0 * (x0\z0); % subtract out projection onto fixed factors
    z2 = bsxfun(@minus, z1, mean(z1));
    [u, s, v] = svds(double(gather(z1)), nF);
    xStart = zscore(u);
end

%w =  zeros(nF+nF0, nDims);

if nargin<6 || isempty(xPrior)
    xPriorConc = [];
else 
    xPriorConc = xPrior*speye(nF);
end

if nargin<7 || isempty(wPrior)
    wPrior = 1e-4;
end

wPriorConc = wPrior*speye(nF+nF0);
wPriorConc(nF+1:nF+nF0,nF+1:nF+nF0)=0;

OldL = -inf;

x = xStart;
w = zeros(nF + nF0, nDims);

for i=1:MaxIter
    % compute weights 
    
    xAll = [x, x0];
    
    [w1, ~, L] = NBreg(y, xAll, r, [], wPriorConc, w);
    if ~isfinite(L)
        fprintf('Infinte score!\n');
        keyboard
    end
    w = w1;

    if Verbose
        fprintf('Iteration %d, score %f\n', i, L);
%         figure(234789502)
%         plot(x(:,1), x(:,2), '.');
%         drawnow
    end

    
    if (L-OldL)<Tol; 
        if (L<OldL); fprintf('NBpca: score went down!\n'); end
        break; 
    end
    OldL = L;
    
    % now recompute scores (not including fixed columns; bias is
    % passed as an offset. remember warm start!
    Offset = (x0 * w(nF+1:nF+nF0,:)) ;
    [x2, ~, L2, MargLike] = NBreg(y', w(1:nF, :)', r, Offset', xPriorConc, x');
    x2 = bsxfun(@minus, x2, mean(x2,2));

    if isempty(xPriorConc)
        % if no x prior, rotate the factors so that x is normalized, w has
        % largest column first
        [u2, s2, v2] = svd(x2', 'econ'); % not including fixed factors
        w2 = s2*v2'*w(1:nF,:)/sqrt(nObs); % and x2 = u2*sqrt(nObs), all x's are orthogonal, meansq 1
        [u3, s3, v3] = svd(w2, 'econ'); 
        w3 = s3*v3';
        x3 = u2*u3*sqrt(nObs);
        % now add back fixed factors and their weights
        x = x3;
        w = [w3; w(nF+1:nF+nF0,:)];
    else
        x = x2'; % just add back fixed factors
    end
        
%     if isempty(xPriorConc)
%         % final iteration to compute marginal likelihood, using a unit
%         % covariance
%         [~, ~, ~, MargLike] = NBreg(y', w(1:nF, :)', r, Offset, speye(nF), x(:,1:nF)');
%     end
%     
%     [x1, ~, L1] = NBreg(y', w', r, [], [], x');
%     x = x1';
    % x = bsxfun(@rdivide, x1', std(x1'));
end

return
%% to test:
N = 1e3;
D = 100;
r = 2;
x0 = [ones(N/2,1), zeros(N/2,1) ; zeros(N/2,1) , ones(N/2,1)];
xReal = [rand(N,1), rand(N,1)]; % 
wReal = [ones(1,D), zeros(1,D) ; zeros(1,D), ones(1,D) ; ...
    3*ones(1,D) , 0*ones(1,D); 0*ones(1,D) , 3*ones(1,D)];
z = [xReal, x0]*wReal;
y = nbinrnd(r, 1 - 1./(1+exp(-z))); % matlab uses other convention for p
[w1, x1] = NBpca(y, r, 2, x0, [], [], 1e-4);%, x);%x);
figure(3987601); clf; plot(x1(:,1), x1(:,2), '.'); % should be a square
%% to test on gene data
r=2;
MyGenes = {'Pvalb', 'Kcnc1', 'Kcna1', 'Slc24a2', 'Gabra1', 'Gabrd', 'Zcchc12', 'Bex1', 'Scg2', '6330403K07Rik', 'Fxyd6', 'Tmsb10', 'Rbp4', 'Tac1', 'Npy'};
MyGenes = mBest.GeneName(mBest.Active);
%gTest = g.CellSubset('.Sst.Pvalb.Rbp4').GeneSubset(MyGenes);
%gTest = g.CellSubset({'.Sst').GoodExpressers(10,10);
%gTest = h.CellSubset({'Pvalb.Tac1'}).ScaleCell(1).GoodExpressers(5,5);
%MyCells = {'Pvalb.Tac1.Syt2', 'Pvalb.Tac1.Nr4a2'};
gTest = h.GoodExpressers(10,5);


y = double(gTest.ScaleCell(1).GeneExp');
[gpIdx, gpName] = grp2idx(gTest.Class);
x0 = dummyvar(gpIdx);
[wCond, xCond] = NBpca(y, r, 1, x0, [], 100, 100);

x0All = ones(gTest.nCells, 1);
[wAll, xAll] = NBpca(y, r, 1, x0, [], 100, 100);
PredMuAll = r*exp([xAll x0All]*wAll);

[gpIdx, gpName] = grp2idx(gTest.TruncateHierarchy(1).Class);
x0Cond1 = dummyvar(gpIdx);
[wCond1, xCond1] = NBpca(y, r, 1, x0Cond1, [], 100, 100);
PredMuCond1 = r*exp([xCond1 x0Cond1]*wCond1);

%%
CGEClass = (1:g.nCells)'>=min(strmatch('Cacna2d1.Ndnf', gTest.Class));
CGERatio = 2*gTest.GeneExp*CGEClass ./ sum(gTest.GeneExp,2) - double(1);

figure(734905); clf; hold on
yCoord = max(wAll(2:end,:), [], 1);
plot(wAll(1,:), yCoord, 'b.');
hand = text(wAll(1,:), yCoord, gTest.GeneName); set(hand, 'color', 'b');
grid on
title('Marginal')
drawnow

figure(734906); clf; hold on
yCoord = max(wCond(2:end,:), [],1);
plot(wCond(1,:), yCoord, 'b.');
hand = text(wCond(1,:), yCoord, gTest.GeneName); set(hand, 'color', 'b');
grid on
title('Class conditional')
drawnow

figure(734907); clf; hold on
yCoord = max(wCond1(2:end,:), [],1);
plot(wCond1(1,:), yCoord, 'b.');
hand = text(wCond1(1,:), yCoord, gTest.GeneName); set(hand, 'color', 'b');
grid on
title('Truncated class conditional')
drawnow
%%
if 1
    MyCells = {'Cck.Lmo1.Vip.Crh', 'Cck.Lmo1.Vip.Fam19a2'};
    figure(734907); 
else
    MyCells = {'Cck.Cxcl14.Slc17a8', 'Cck.Cxcl14.Calb1.Kctd12', 'Cck.Cxcl14.Calb1.Npy', 'Cck.Cxcl14.Calb1.Tnfaip8l3'};
    figure(734908); 
end

gClass = h.CellSubset(MyCells).GoodExpressers(5,2);
CellIDs = h.IdentifyCells(MyCells);
[wClass, xClass] = NBpca(double(gClass.ScaleCell(1).GeneExp'), r, 1, [], [], .1, .1);
%[wClass, xClass] = NBpca(double(gClass.ScaleCell(1).GeneExp'), r, 1, [], xAll(CellIDs), 1, 1);
clf; hold on
yCoord = sum(wClass(2:end,:), 1);
plot(wClass(1,:), yCoord, 'b.');
hand = text(wClass(1,:), yCoord, gClass.GeneName); set(hand, 'color', 'b');
grid on
title(sprintf('%s ', MyCells{:}));
drawnow

% PredMu = r*exp([x x0]*w);
%%
MyGeneName = 'Rgs10';
gi = gTest.NamesToIDs(MyGeneName);
figure(98704365)
%gTest.ScaleCell(1).Scatter(PredMuCond1(:,gi), gi);
gTest.ScaleCell(1).Scatter(gi, PredMuCond1(:,gi));
figure(98704366)
gTest.ScaleCell(1).Scatter(gi, PredMuAll(:,gi));


%%

