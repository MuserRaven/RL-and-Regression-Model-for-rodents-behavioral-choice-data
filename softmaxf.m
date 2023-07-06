function y=softmaxf(A,B,beta,bias)
% A=qR, B=qL
%y=1./(1+exp(-(bias+beta)*(A-B))); 
 y=exp((A+bias)*beta)./(exp(B*beta)+ exp((A+bias)*beta));  % the probability that you choose A 
% for all A,B,beta, and bias belongs to real number 0<y<1

% compared to the ideal, if bias is postive, it means you are underweighting the proportion of A,
% if negative, it means you are overweighting the proportion of A 


end

