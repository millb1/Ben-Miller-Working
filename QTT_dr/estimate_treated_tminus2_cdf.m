function [fval] =estimate_treated_tminus2_cdf(Y1_treat1,delta1)


indicator_Y1treat =double((Y1_treat1 <= delta1));

fval = mean(indicator_Y1treat,1);
