function [wAll, pAll, SumL, MargLike, Hout] = NBreg(yAll, x, r, OffsetAll, PriorConc, w0);
% [w, p, L, MargLike, Hess] = NBreg(y, x, r, Offset, PriorConc);
% 
% does negative binomial regression to predict y from x. Each row of y and
% x represents an observation - each column of y is predicted
% independently. Don't forget to include a constant column in x (with 
% similar norm to other comlumns, for numerical reasons).
% Offset is an optional matrix added to x*w (broadcast if necessary)
% PriorConc is for the weight vectors (assuming mean 0 for each one)
% w0 is a warm start
%
% y is size nObs by nVars
% x is size nObs by nPreds
% output w is nPreds by nVar weight matrix
% optional output p is nObs by nVars set of predicted probabilities
% optional output L is log likelihood
% optional output MarkLike gives the marginal likelihood for each variable,
% optional output Hess is the Hessian matrix of the likielihood at the optimal value 
%
% MODEL:
% p = 1./(1+exp(-(x*w + Offset)))
% L = const + r*log(1-p) + y.*log(p)
% So predicted mean yhat = r*exp(x*w)

[nObs, nVar] = size(yAll);
[nObs2, nPred] = size(x);

if nargin<4 || isempty(OffsetAll) ; OffsetAll=0; end;
if nargin<5 || isempty(PriorConc); PriorConc = zeros(nPred); end;
if nargin<6 || isempty(w0); w0 = zeros(nPred, nVar); end

MaxIter = 100;
Tol = 1e-3;
Verbose = 0;
Reg = 0; % add before inverting Hessian.


if nObs~=nObs2
    error('x and y should have same number of rows');
end

wAll = zeros(nPred, nVar);
pAll = zeros(nObs, nVar);
MargLike = zeros(nVar,1);

SumL = 0;
if nargout>=5
    Hout = zeros(nPred,nPred,nVar);
end
for v=1:nVar
    nShrinkSteps = 0;
    y = yAll(:,v);
    w = w0(:,v);
    if size(OffsetAll,2)>1
        Offset = OffsetAll(:,v); 
    else
        Offset = OffsetAll;
    end
    OldL = -inf;
    ExpectedIncrease = 0;
    for i=1:MaxIter
        p = 1./(1+exp(-bsxfun(@plus, x*w, Offset))); % nObs by 1
        L = r*sum(log(1-p(:))) + y'*log(p); % scalar likelihood
        
        if L<OldL+ExpectedIncrease
            % score went down - so cut the step by half and try again
            dw = dw/2;
            w = w-dw;
            ExpectedIncrease = ExpectedIncrease/2;
            nShrinkSteps = nShrinkSteps+1;
            continue;
        end
        
        if Verbose==2
            fprintf('NBreg: Var %d, Iter %d, likelihood %.4f\n', v, i, L);
        end
        
        if (L-OldL)<Tol
            if L<OldL
                fprintf('NBreg: score went down\n');
            end
            break;
        end
        OldL = L;
        
%         if nShrinkSteps>0
%             fprintf('NBreg: continuing after shrink\n');
%         end
        
        % compute gradient
        dLdw = (y - p.*(y+r))'*x - w'*PriorConc; 
        % negative Hessian
        H = x' * bsxfun(@times, x, (y+r).*p.*(1-p)) + PriorConc;
        % update
        dw = (H + Reg*eye(nPred)) \ dLdw';
        w = w + dw;
        ExpectedIncrease = dLdw*dw*.1;
    end
    if Verbose && (i==MaxIter || nShrinkSteps>0)
        fprintf('NBreg: %d iterations, %d shrink steps\n', i, nShrinkSteps);
    end
        
        
    wAll(:,v) = gather(w);
    pAll(:,v) = gather(p);
    SumL = SumL+L;
    % marginal likelihood: Bishop p. 216
    if nargout>=4
        MargLike(v) = L - w'*PriorConc*w/2 + log(det(PriorConc))/2 - log(det(H))/2;
    end
    if nargout>=5
        Hout(:,:,v) = H;
    end
end
        
return
%% to test:
N = 1e4;
r = 5;
x = [randn(N,2), ones(N,1)]; % each column has rms of 1
wReal = [1e0*[1 0 ; 0 1]; 5 5];
z = x*wReal;
y = nbinrnd(r, 1 - 1./(1+exp(-z))); % matlab uses other convention for p
[w p L Hess] = NBreg(y, x, r, ones(N,1) * [2 3]);
w
return
%% now estimate points from weights:
[xEst p2 L2 MargLike] = NBreg(y',wReal(1:2,:)', r, wReal(3,:)', eye(2));
%[xEst3 p3 L3 Hess3] = NBreg(y',wReal', r);
figure(29034857); clf; hold on
plot(x(:,1), xEst(1,:), 'b.')
plot(x(:,2), xEst(2,:), 'r.')
grid on;
plot(xlim, xlim, 'k:');

%% and compute likelihood correction
L1 = zeros(N,1);
L0 = zeros(N,1);
for i=1:N
    H = Hess3(:,:,i);
    L1(i) = - .5*log(det(eye(2) + H)) - .5*xEst(:,i)'*inv(eye(2)+inv(H))*xEst(:,i);
    L0(i) = - .5*xEst(:,i)'*xEst(:,i);
end

%% To test on real data:
gTest = g.CellSubset('.Sst.Pvalb.Pvalb').GoodExpressers(4,30);
[~, gPred] = gTest.Log(30).PCA(2);
x = [zscore(gPred.GeneExp'), ones(gPred.nCells, 1)];
%x = [ones(gPred.nCells, 1)];
y = gTest.GeneExp';
w = NBreg(y, x, 2);

PredMu = 2*exp(x*w);

MyNames = {'Pvalb', 'Kcna1', 'Impact'};

figure(7001);

eg = gTest.NamesToIDs(MyNames);
plot3(y(:,eg(1)), y(:,eg(2)), y(:,eg(3)), 'b.', PredMu(:,eg(1)), PredMu(:,eg(2)), PredMu(:,eg(3)), 'r.');
xlabel(MyNames{1});
ylabel(MyNames{2});
zlabel(MyNames{3});

% now given the weights, find where cells live on a continuum
[xEst p2 L2 MargLike] = NBreg(y',w(1:2,:)', r, w(3,:)', eye(2));

