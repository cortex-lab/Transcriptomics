function Lse = LogSumExp(L);
% p = LogLtoP(L)
%
% given a set of log likelihoods L, robustly compute the log sum exp:
% L1(i,c) = log( sum_i(exp(L(i,c)));
% operates column-wise

MaxL = max(L,[],1);
L1 = bsxfun(@minus,L, MaxL);
eL = exp(L1);
Lse = bsxfun(@plus, log(sum(exp(L1),1)), MaxL);