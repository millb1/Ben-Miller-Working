function [fval] =estimate_treated_dtminus1_cdf(DY11,delta2)


indicator_DY1treat =double((DY11 <= delta2));

fval = mean(indicator_DY1treat,1);
