function [fval] =estimate_untreateddiff_treated_alt_cdf_procorrect(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,prob_weight,delta)

L=exp(gammahat(1,1)*x1 + gammahat(1,2)*x2 +gammahat(1,3)*x3+ gammahat(1,4)*x4)./(1+exp(gammahat(1,1)*x1 + gammahat(1,2)*x2 +gammahat(1,3)*x3 + gammahat(1,4)*x4));

f_element = zeros(size(DY2));

indicator_DY2 =double((DY2 <= delta));

d=betahat(1,1)*x0 + betahat(2,1)*x1 + betahat(3,1)*x2 + betahat(4,1)*x3 + betahat(5,1)*x4;

for p=1:size(DY2)
    f_element(p,1) =((((1-D(p,1))/prob_weight)*(L(p,1)/(1-L(p,1))))*indicator_DY2(p,1))-((((1-D(p,1))/prob_weight)*(L(p,1)/(1-L(p,1)))-(D(p,1)/prob))*unifcdf(delta - d(p,1)/2));
end

fval = mean(f_element,1);