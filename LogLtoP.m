function p = LogLtoP(L);
% p = LogLtoP(L)
%
% given a set of log likelihoods L, robustly compute the posterior
% probabilities p(i,c) = exp(L(i,c))/ sum_i(exp(L(i,c)));
% operates column-wise

L1 = bsxfun(@minus,L, max(L,[],1));
eL = exp(L1);
p = bsxfun(@rdivide, eL, sum(eL,1));