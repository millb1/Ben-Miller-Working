function[fval] = estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input)

indicator_Y0treat=zeros(size(q_tminus1));
for p=1:size(q_tminus1)
    indicator_Y0treat(p,1) = double(q_tminus1_0treated(p,1) <= y_input -q_tminus1(p,1));
end
fval = mean(indicator_Y0treat,1);
