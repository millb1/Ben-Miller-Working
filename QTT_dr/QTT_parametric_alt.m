N=5000;
T=3;
nloop = 100;
rng(1);
tau = 0.5;
x0_full = ones(N,nloop);
x1_full = normrnd(0,1,[N,nloop]);
x2_full = normrnd(0,2,[N,nloop]);
x3_full = normrnd(0,3,[N,nloop]);
x4_full = normrnd(0,4,[N,nloop]);

iqtt = zeros(nloop,1);
%both first stage estimates are correct
for k=1:nloop
fprintf('Iterations No %i ',k);
x0 = x0_full(:,k);
x1 =x1_full(:,k);
x2 =x2_full(:,k);
x3 =x3_full(:,k);
x4 =x4_full(:,k);

pro = exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4)./(1+exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4));

D=zeros(size(x1));
for p=1:size(D)
D(p,1)=binornd(1,pro(p,1));
end
eta = zeros(size(D));
for p = 1:size(D)
    if D(p,1) ==0
eta(p,1)= normrnd(0,1); 
    else
    eta(p,1)= normrnd(1,1);
    end
end
v = normrnd(0,1,[N,T]);
Y1 = 0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,1);
Y2 =1 +0.50*x1 + 0.75*x2 +1.00*x3 + 1.5*x4 + eta + v(:,2);
Y3_1 =2*(0.75*x1 + 0.5*x2 +0.75*x3 + 1*x4) + eta + v(:,3);
Y3_0 =(2 +0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,3));
Y3= D.*Y3_1 +(1-D).*Y3_0;
DY1 = Y2 -Y1;
DY2 = Y3-Y2;

Y3_0_1 = zeros(size(D));
for p=1:size(Y3_0_1)
    if D(p,1)==0
        Y3_0_1 (p,1) = 0;
    else
       Y3_0_1 (p,1) = Y3_0(p,1); 
    end
end
Y3_01 = Y3_0_1(Y3_0_1~=0);

DY1_1 = zeros(size(D));
for p=1:size(DY1_1)
    if D(p,1) ==0
        DY1_1(p,1) =0;
    else
        DY1_1(p,1) = DY1(p,1);
    end
end
DY11 = DY1_1(DY1_1~=0);

DY2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_0(p,1) =DY2(p,1);
    else
        DY2_0(p,1) = 0;
    end
end
DY20 = DY2_0(DY2_0~=0);

DY2_1 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_1(p,1) =0;
    else
        DY2_1(p,1) =DY2(p,1);
    end
end
DY21 = DY2_1(DY2_1~=0);

x0_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x0_0(p,1) =x0(p,1);
    else
        x0_0(p,1) = 0;
    end
end
x00 = x0_0(x0_0~=0);

x1_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x1_0(p,1) =x1(p,1);
    else
        x1_0(p,1) = 0;
    end
end
x10 = x1_0(x1_0~=0);

x2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x2_0(p,1) =x2(p,1);
    else
        x2_0(p,1) = 0;
    end
end
x20 = x2_0(x2_0~=0);

x3_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x3_0(p,1) =x3(p,1);
    else
        x3_0(p,1) = 0;
    end
end
x30 = x3_0(x3_0~=0);

x4_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x4_0(p,1) =x4(p,1);
    else
        x4_0(p,1) = 0;
    end
end
x40 = x4_0(x4_0~=0);

Y1_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y1_treat(p,1) =0;
    else
        Y1_treat(p,1) =Y1(p,1);
    end
end
Y1_treat1=Y1_treat(Y1_treat~=0);

Y2_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y2_treat(p,1) =0;
    else
        Y2_treat(p,1) =Y2(p,1);
    end
end
Y2_treat1=Y2_treat(Y2_treat~=0);

Y3_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y3_treat(p,1) =0;
    else
        Y3_treat(p,1) =Y3(p,1);
    end
end
Y3_treat1=Y3_treat(Y3_treat~=0);





%First Stage: Estimate beta for conditional cdf

X_data = [x00(:),x10(:),x20(:),x30(:),x40(:)];

betahat = (X_data'*X_data)\(X_data'*DY20);

%First Stage: Estimate logit
%initial = [0,-0.25,-0.5,-0.75,-1]; %Estimate gammas
%W = eye(5);
%options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
%[gammahat] = fminsearch(@(gamma) GMM_pro(gamma,x0,x1,x2,x3,x4,D,W),initial,options);
mynegloglik = @(gamma) -sum(log(binopdf(D,1,exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)./(1+exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)))));
initial = [-0.25,-0.5,-0.75,-1];
opts = optimset('fminsearch');
opts.MaxFunEvals = Inf;
opts.MaxIter = 10000;
gammahat = fminsearch(mynegloglik,initial,opts);


probx =exp(gammahat(1)*x1 + gammahat(2)*x2 + gammahat(3)*x3 +gammahat(4)*x4)./(1+exp(gammahat(1)*x1 + gammahat(2)*x2 + gammahat(3)*x3 +gammahat(4)*x4)); 

weight = ((1-D).*probx)./(1-probx);

prob_weight = mean(weight,1);



%Generate empirical cdfs and cdf estimates
[f,x_y3treat] = ecdf(Y3_treat1);
q_1 = interp1(f,x_y3treat,tau,'next');

prob =mean(D);
%delta=0 This is a value to test the function below. Also, generating
%doubly robust portion
%estimate_untreateddiff_treated=estimate_untreateddiff_treated_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,delta);



%estimate_tminus2_treated = estimate_treated_tminus2_cdf(Y1_treat1,delta1);

%estimate_dtminus1_treated = estimate_treated_dtminus1_cdf(DY21,delta2);

estimate_tminus2_treated =zeros(size(Y1_treat1));
%generating empirical cdf at values of Y1_treat1 (using the first period data of the treated) 
for p=1:size(Y1_treat1)
    indicator_Y1treat =double((Y1_treat1 <= Y1_treat1(p,1)));
    estimate_tminus2_treated(p,1)= mean(indicator_Y1treat,1);
end

[f1,x_y2treat] = ecdf(Y2_treat1);
q_tminus1 = zeros(size(Y2_treat1));
for p=1:size(Y2_treat1)
    q_tminus1(p,1) = interp1(f1,x_y2treat,estimate_tminus2_treated(p,1),'next');
end

estimate_dtminus1_treated = zeros(size(DY11));
for p=1:size(DY11)
    indicator_DY11 = double((DY11 <= DY11(p,1)));
    estimate_dtminus1_treated(p,1)=mean(indicator_DY11,1);
end

x_deltaY0t = -12:.001:12;
x_deltaY0t =x_deltaY0t';

%estimate_untreateddiff_treated=zeros(size(estimate_dtminus1_treated));
q_tminus1_0treated = zeros(size(estimate_dtminus1_treated));
for p=1:size(estimate_dtminus1_treated)
 estimate_untreateddiff_treated1=@(delta) estimate_untreateddiff_treated_alt_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,prob_weight,delta)-estimate_dtminus1_treated(p,1);
q_tminus1_0treated(p,1) = fzero(estimate_untreateddiff_treated1,0);
end
estimate_untreated_treated_q = @(y_input) estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input) - tau;
q_0 = fzero(estimate_untreated_treated_q,0);


iqtt(k,1) = q_1-q_0;
end
qtt_mean = mean(iqtt,1);

qtt_std = std(iqtt);

%propensity score correct
N=1000;
T=3;
nloop = 200;
rng(1);
tau = 0.75;
x0_full = ones(N,nloop);
x1_full = normrnd(0,1,[N,nloop]);
x2_full = normrnd(0,2,[N,nloop]);
x3_full = normrnd(0,3,[N,nloop]);
x4_full = normrnd(0,4,[N,nloop]);


for k=1:nloop
fprintf('Iterations No %i ',k);
x0 = x0_full(:,k);
x1 =x1_full(:,k);
x2 =x2_full(:,k);
x3 =x3_full(:,k);
x4 =x4_full(:,k);

pro = exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4)./(1+exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4));

D=zeros(size(x1));
for p=1:size(D)
D(p,1)=binornd(1,pro(p,1));
end
eta = zeros(size(D));
for p = 1:size(D)
    if D(p,1) ==0
eta(p,1)= normrnd(0,1); 
    else
    eta(p,1)= normrnd(1,1);
    end
end
v = normrnd(0,1,[N,T]);
Y1 = 0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,1);
Y2 =1 +0.50*x1 + 0.75*x2 +1.00*x3 + 1.5*x4 + eta + v(:,2);
Y3_1 =2*(0.75*x1 + 0.5*x2 +0.75*x3 + 1*x4) + eta + v(:,3);
Y3_0 =(2 +0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,3));
Y3= D.*Y3_1 +(1-D).*Y3_0;
DY1 = Y2 -Y1;
DY2 = Y3-Y2;

Y3_0_1 = zeros(size(D));
for p=1:size(Y3_0_1)
    if D(p,1)==0
        Y3_0_1 (p,1) = 0;
    else
       Y3_0_1 (p,1) = Y3_0(p,1); 
    end
end
Y3_01 = Y3_0_1(Y3_0_1~=0);

DY1_1 = zeros(size(D));
for p=1:size(DY1_1)
    if D(p,1) ==0
        DY1_1(p,1) =0;
    else
        DY1_1(p,1) = DY1(p,1);
    end
end
DY11 = DY1_1(DY1_1~=0);

DY2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_0(p,1) =DY2(p,1);
    else
        DY2_0(p,1) = 0;
    end
end
DY20 = DY2_0(DY2_0~=0);

DY2_1 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_1(p,1) =0;
    else
        DY2_1(p,1) =DY2(p,1);
    end
end
DY21 = DY2_1(DY2_1~=0);

x0_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x0_0(p,1) =x0(p,1);
    else
        x0_0(p,1) = 0;
    end
end
x00 = x0_0(x0_0~=0);

x1_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x1_0(p,1) =x1(p,1);
    else
        x1_0(p,1) = 0;
    end
end
x10 = x1_0(x1_0~=0);

x2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x2_0(p,1) =x2(p,1);
    else
        x2_0(p,1) = 0;
    end
end
x20 = x2_0(x2_0~=0);

x3_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x3_0(p,1) =x3(p,1);
    else
        x3_0(p,1) = 0;
    end
end
x30 = x3_0(x3_0~=0);

x4_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x4_0(p,1) =x4(p,1);
    else
        x4_0(p,1) = 0;
    end
end
x40 = x4_0(x4_0~=0);

Y1_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y1_treat(p,1) =0;
    else
        Y1_treat(p,1) =Y1(p,1);
    end
end
Y1_treat1=Y1_treat(Y1_treat~=0);

Y2_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y2_treat(p,1) =0;
    else
        Y2_treat(p,1) =Y2(p,1);
    end
end
Y2_treat1=Y2_treat(Y2_treat~=0);

Y3_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y3_treat(p,1) =0;
    else
        Y3_treat(p,1) =Y3(p,1);
    end
end
Y3_treat1=Y3_treat(Y3_treat~=0);





%First Stage: Estimate beta for conditional cdf

X_data = [x00(:),x10(:),x20(:),x30(:),x40(:)];

betahat = (X_data'*X_data)\(X_data'*DY20);

%First Stage: Estimate logit
%initial = [0,-0.25,-0.5,-0.75,-1]; %Estimate gammas
%W = eye(5);
%options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
%[gammahat] = fminsearch(@(gamma) GMM_pro(gamma,x0,x1,x2,x3,x4,D,W),initial,options);
mynegloglik = @(gamma) -sum(log(binopdf(D,1,exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)./(1+exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)))));
initial = [-0.25,-0.5,-0.75,-1];
opts = optimset('fminsearch');
opts.MaxFunEvals = Inf;
opts.MaxIter = 10000;
gammahat = fminsearch(mynegloglik,initial,opts);

probx =exp(gammahat(1)*x1 + gammahat(2)*x2 + gammahat(3)*x3 +gammahat(4)*x4)./(1+exp(gammahat(1)*x1 + gammahat(2)*x2 + gammahat(3)*x3 +gammahat(4)*x4)); 

weight = ((1-D).*probx)./(1-probx);

prob_weight = mean(weight,1);

%Generate empirical cdfs and cdf estimates
[f,x_y3treat] = ecdf(Y3_treat1);
q_1 = interp1(f,x_y3treat,tau,'next');

prob =mean(D);
%delta=0 This is a value to test the function below. Also, generating
%doubly robust portion
%estimate_untreateddiff_treated=estimate_untreateddiff_treated_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,delta);



%estimate_tminus2_treated = estimate_treated_tminus2_cdf(Y1_treat1,delta1);

%estimate_dtminus1_treated = estimate_treated_dtminus1_cdf(DY21,delta2);

estimate_tminus2_treated =zeros(size(Y1_treat1));
%generating empirical cdf at values of Y1_treat1 (using the first period data of the treated) 
for p=1:size(Y1_treat1)
    indicator_Y1treat =double((Y1_treat1 <= Y1_treat1(p,1)));
    estimate_tminus2_treated(p,1)= mean(indicator_Y1treat,1);
end

[f1,x_y2treat] = ecdf(Y2_treat1);
q_tminus1 = zeros(size(Y2_treat1));
for p=1:size(Y2_treat1)
    q_tminus1(p,1) = interp1(f1,x_y2treat,estimate_tminus2_treated(p,1),'next');
end

estimate_dtminus1_treated = zeros(size(DY11));
for p=1:size(DY11)
    indicator_DY11 = double((DY11 <= DY11(p,1)));
    estimate_dtminus1_treated(p,1)=mean(indicator_DY11,1);
end

x_deltaY0t = -12:.001:12;
x_deltaY0t =x_deltaY0t';

%estimate_untreateddiff_treated=zeros(size(estimate_dtminus1_treated));
q_tminus1_0treated = zeros(size(estimate_dtminus1_treated));
for p=1:size(estimate_dtminus1_treated)
 estimate_untreateddiff_treated1=@(delta) estimate_untreateddiff_treated_alt_cdf_procorrect(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,prob_weight,delta)-estimate_dtminus1_treated(p,1);
q_tminus1_0treated(p,1) = fzero(estimate_untreateddiff_treated1,0);
end
estimate_untreated_treated_q = @(y_input) estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input) - tau;
q_0 = fzero(estimate_untreated_treated_q,0);


iqtt(k,1) = q_1-q_0;
end
qtt_mean = mean(iqtt,1);

qtt_std = std(iqtt);


%Conditional CDF correct


N=1000;
T=3;
nloop = 200;
rng(1);
tau = 0.75;
x0_full = ones(N,nloop);
x1_full = normrnd(0,1,[N,nloop]);
x2_full = normrnd(0,2,[N,nloop]);
x3_full = normrnd(0,3,[N,nloop]);
x4_full = normrnd(0,4,[N,nloop]);

iqtt = zeros(nloop,1);
for k=1:nloop
fprintf('Iterations No %i ',k);
x0 = x0_full(:,k);
x1 =x1_full(:,k);
x2 =x2_full(:,k);
x3 =x3_full(:,k);
x4 =x4_full(:,k);

pro = exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4)./(1+exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4));

D=zeros(size(x1));
for p=1:size(D)
D(p,1)=binornd(1,pro(p,1));
end
eta = zeros(size(D));
for p = 1:size(D)
    if D(p,1) ==0
eta(p,1)= normrnd(0,1); 
    else
    eta(p,1)= normrnd(1,1);
    end
end
v = normrnd(0,1,[N,T]);
Y1 = 0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,1);
Y2 =1 +0.50*x1 + 0.75*x2 +1.00*x3 + 1.5*x4 + eta + v(:,2);
Y3_1 =2*(0.75*x1 + 0.5*x2 +0.75*x3 + 1*x4) + eta + v(:,3);
Y3_0 =(2 +0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,3));
Y3= D.*Y3_1 +(1-D).*Y3_0;
DY1 = Y2 -Y1;
DY2 = Y3-Y2;

Y3_0_1 = zeros(size(D));
for p=1:size(Y3_0_1)
    if D(p,1)==0
        Y3_0_1 (p,1) = 0;
    else
       Y3_0_1 (p,1) = Y3_0(p,1); 
    end
end
Y3_01 = Y3_0_1(Y3_0_1~=0);

DY1_1 = zeros(size(D));
for p=1:size(DY1_1)
    if D(p,1) ==0
        DY1_1(p,1) =0;
    else
        DY1_1(p,1) = DY1(p,1);
    end
end
DY11 = DY1_1(DY1_1~=0);

DY2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_0(p,1) =DY2(p,1);
    else
        DY2_0(p,1) = 0;
    end
end
DY20 = DY2_0(DY2_0~=0);

DY2_1 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_1(p,1) =0;
    else
        DY2_1(p,1) =DY2(p,1);
    end
end
DY21 = DY2_1(DY2_1~=0);

x0_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x0_0(p,1) =x0(p,1);
    else
        x0_0(p,1) = 0;
    end
end
x00 = x0_0(x0_0~=0);

x1_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x1_0(p,1) =x1(p,1);
    else
        x1_0(p,1) = 0;
    end
end
x10 = x1_0(x1_0~=0);

x2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x2_0(p,1) =x2(p,1);
    else
        x2_0(p,1) = 0;
    end
end
x20 = x2_0(x2_0~=0);

x3_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x3_0(p,1) =x3(p,1);
    else
        x3_0(p,1) = 0;
    end
end
x30 = x3_0(x3_0~=0);

x4_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x4_0(p,1) =x4(p,1);
    else
        x4_0(p,1) = 0;
    end
end
x40 = x4_0(x4_0~=0);

Y1_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y1_treat(p,1) =0;
    else
        Y1_treat(p,1) =Y1(p,1);
    end
end
Y1_treat1=Y1_treat(Y1_treat~=0);

Y2_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y2_treat(p,1) =0;
    else
        Y2_treat(p,1) =Y2(p,1);
    end
end
Y2_treat1=Y2_treat(Y2_treat~=0);

Y3_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y3_treat(p,1) =0;
    else
        Y3_treat(p,1) =Y3(p,1);
    end
end
Y3_treat1=Y3_treat(Y3_treat~=0);





%First Stage: Estimate beta for conditional cdf

X_data = [x00(:),x10(:),x20(:),x30(:),x40(:)];

betahat = (X_data'*X_data)\(X_data'*DY20);

%First Stage: Estimate logit
%initial = [0,-0.25,-0.5,-0.75,-1]; %Estimate gammas
%W = eye(5);
%options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
%[gammahat] = fminsearch(@(gamma) GMM_pro(gamma,x0,x1,x2,x3,x4,D,W),initial,options);
%mynegloglik = @(gamma) -sum(log(binopdf(D,1,cdf('chi2',gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4,1))));
%mynegloglik = @(gamma) -sum(log(binopdf(D,1,chi2cdf(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4,1))));
%initial = [-0.25,-0.5,-0.75,-1];
%opts = optimset('fminsearch');
%opts.MaxFunEvals = Inf;
%opts.MaxIter = 10000;
%gammahat = fminsearch(mynegloglik,initial,opts);

%probx =cdf('chi2',gammahat(1)*x1 + gammahat(2)*x2 + gammahat(3)*x3 +gammahat(4)*x4); 
mynegloglik = @(gamma) -sum(log(binopdf(D,1,exp(gamma(1)*x1)./(1+exp(gamma(1)*x1)))));
initial = [-0.25];
opts = optimset('fminsearch');
opts.MaxFunEvals = Inf;
opts.MaxIter = 10000;
gammahat = fminsearch(mynegloglik,initial,opts);

probx =exp(gammahat(1)*x1)./(1+exp(gammahat(1)*x1)); 



weight = ((1-D).*probx)./(1-probx);

prob_weight = mean(weight,1);


%Generate empirical cdfs and cdf estimates
[f,x_y3treat] = ecdf(Y3_treat1);
q_1 = interp1(f,x_y3treat,tau,'next');

prob =mean(D);
%delta=0 This is a value to test the function below. Also, generating
%doubly robust portion
%estimate_untreateddiff_treated=estimate_untreateddiff_treated_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,delta);



%estimate_tminus2_treated = estimate_treated_tminus2_cdf(Y1_treat1,delta1);

%estimate_dtminus1_treated = estimate_treated_dtminus1_cdf(DY21,delta2);

estimate_tminus2_treated =zeros(size(Y1_treat1));
%generating empirical cdf at values of Y1_treat1 (using the first period data of the treated) 
for p=1:size(Y1_treat1)
    indicator_Y1treat =double((Y1_treat1 <= Y1_treat1(p,1)));
    estimate_tminus2_treated(p,1)= mean(indicator_Y1treat,1);
end

[f1,x_y2treat] = ecdf(Y2_treat1);
q_tminus1 = zeros(size(Y2_treat1));
for p=1:size(Y2_treat1)
    q_tminus1(p,1) = interp1(f1,x_y2treat,estimate_tminus2_treated(p,1),'next');
end

estimate_dtminus1_treated = zeros(size(DY11));
for p=1:size(DY11)
    indicator_DY11 = double((DY11 <= DY11(p,1)));
    estimate_dtminus1_treated(p,1)=mean(indicator_DY11,1);
end

x_deltaY0t = -12:.001:12;
x_deltaY0t =x_deltaY0t';

%estimate_untreateddiff_treated=zeros(size(estimate_dtminus1_treated));
q_tminus1_0treated = zeros(size(estimate_dtminus1_treated));
for p=1:size(estimate_dtminus1_treated)
 estimate_untreateddiff_treated1=@(delta) estimate_untreateddiff_treated_alt_cdf_condcdfcorrect(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,prob_weight,delta)-estimate_dtminus1_treated(p,1);
q_tminus1_0treated(p,1) = fzero(estimate_untreateddiff_treated1,1);
end
estimate_untreated_treated_q = @(y_input) estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input) - tau;
q_0 = fzero(estimate_untreated_treated_q,1);


iqtt(k,1) = q_1-q_0;
end
qtt_mean = mean(iqtt,1);

qtt_std = std(iqtt);




%Callaway and Li Estimator

N=1000;
T=3;
nloop = 200;
rng(1);
tau = 0.5;
x0_full = ones(N,nloop);
x1_full = normrnd(0,1,[N,nloop]);
x2_full = normrnd(0,2,[N,nloop]);
x3_full = normrnd(0,3,[N,nloop]);
x4_full = normrnd(0,4,[N,nloop]);

iqtt = zeros(nloop,1);
%both first stage estimates are correct
for k=1:nloop
fprintf('Iterations No %i ',k);
x0 = x0_full(:,k);
x1 =x1_full(:,k);
x2 =x2_full(:,k);
x3 =x3_full(:,k);
x4 =x4_full(:,k);

pro = exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4)./(1+exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4));

D=zeros(size(x1));
for p=1:size(D)
D(p,1)=binornd(1,pro(p,1));
end
eta = zeros(size(D));
for p = 1:size(D)
    if D(p,1) ==0
eta(p,1)= normrnd(0,1); 
    else
    eta(p,1)= normrnd(1,1);
    end
end
v = normrnd(0,1,[N,T]);
Y1 = 0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,1);
Y2 =1 +0.50*x1 + 0.75*x2 +1.00*x3 + 1.5*x4 + eta + v(:,2);
Y3_1 =2*(0.75*x1 + 0.5*x2 +0.75*x3 + 1*x4) + eta + v(:,3);
Y3_0 =(2 +0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,3));
Y3= D.*Y3_1 +(1-D).*Y3_0;
DY1 = Y2 -Y1;
DY2 = Y3-Y2;

Y3_0_1 = zeros(size(D));
for p=1:size(Y3_0_1)
    if D(p,1)==0
        Y3_0_1 (p,1) = 0;
    else
       Y3_0_1 (p,1) = Y3_0(p,1); 
    end
end
Y3_01 = Y3_0_1(Y3_0_1~=0);

DY1_1 = zeros(size(D));
for p=1:size(DY1_1)
    if D(p,1) ==0
        DY1_1(p,1) =0;
    else
        DY1_1(p,1) = DY1(p,1);
    end
end
DY11 = DY1_1(DY1_1~=0);

DY2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_0(p,1) =DY2(p,1);
    else
        DY2_0(p,1) = 0;
    end
end
DY20 = DY2_0(DY2_0~=0);

DY2_1 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_1(p,1) =0;
    else
        DY2_1(p,1) =DY2(p,1);
    end
end
DY21 = DY2_1(DY2_1~=0);

x0_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x0_0(p,1) =x0(p,1);
    else
        x0_0(p,1) = 0;
    end
end
x00 = x0_0(x0_0~=0);

x1_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x1_0(p,1) =x1(p,1);
    else
        x1_0(p,1) = 0;
    end
end
x10 = x1_0(x1_0~=0);

x2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x2_0(p,1) =x2(p,1);
    else
        x2_0(p,1) = 0;
    end
end
x20 = x2_0(x2_0~=0);

x3_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x3_0(p,1) =x3(p,1);
    else
        x3_0(p,1) = 0;
    end
end
x30 = x3_0(x3_0~=0);

x4_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x4_0(p,1) =x4(p,1);
    else
        x4_0(p,1) = 0;
    end
end
x40 = x4_0(x4_0~=0);

Y1_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y1_treat(p,1) =0;
    else
        Y1_treat(p,1) =Y1(p,1);
    end
end
Y1_treat1=Y1_treat(Y1_treat~=0);

Y2_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y2_treat(p,1) =0;
    else
        Y2_treat(p,1) =Y2(p,1);
    end
end
Y2_treat1=Y2_treat(Y2_treat~=0);

Y3_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y3_treat(p,1) =0;
    else
        Y3_treat(p,1) =Y3(p,1);
    end
end
Y3_treat1=Y3_treat(Y3_treat~=0);





%First Stage: Estimate beta for conditional cdf


%First Stage: Estimate logit
%initial = [0,-0.25,-0.5,-0.75,-1]; %Estimate gammas
%W = eye(5);
%options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
%[gammahat] = fminsearch(@(gamma) GMM_pro(gamma,x0,x1,x2,x3,x4,D,W),initial,options);
mynegloglik = @(gamma) -sum(log(binopdf(D,1,exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)./(1+exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)))));
initial = [-0.25,-0.5,-0.75,-1];
opts = optimset('fminsearch');
opts.MaxFunEvals = Inf;
opts.MaxIter = 10000;
gammahat = fminsearch(mynegloglik,initial,opts);


%Generate empirical cdfs and cdf estimates
[f,x_y3treat] = ecdf(Y3_treat1);
q_1 = interp1(f,x_y3treat,tau,'next');

prob =mean(D);
%delta=0 This is a value to test the function below. Also, generating
%doubly robust portion
%estimate_untreateddiff_treated=estimate_untreateddiff_treated_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,delta);



%estimate_tminus2_treated = estimate_treated_tminus2_cdf(Y1_treat1,delta1);

%estimate_dtminus1_treated = estimate_treated_dtminus1_cdf(DY21,delta2);

estimate_tminus2_treated =zeros(size(Y1_treat1));
%generating empirical cdf at values of Y1_treat1 (using the first period data of the treated) 
for p=1:size(Y1_treat1)
    indicator_Y1treat =double((Y1_treat1 <= Y1_treat1(p,1)));
    estimate_tminus2_treated(p,1)= mean(indicator_Y1treat,1);
end

[f1,x_y2treat] = ecdf(Y2_treat1);
q_tminus1 = zeros(size(Y2_treat1));
for p=1:size(Y2_treat1)
    q_tminus1(p,1) = interp1(f1,x_y2treat,estimate_tminus2_treated(p,1),'next');
end

estimate_dtminus1_treated = zeros(size(DY11));
for p=1:size(DY11)
    indicator_DY11 = double((DY11 <= DY11(p,1)));
    estimate_dtminus1_treated(p,1)=mean(indicator_DY11,1);
end

x_deltaY0t = -12:.001:12;
x_deltaY0t =x_deltaY0t';

%estimate_untreateddiff_treated=zeros(size(estimate_dtminus1_treated));
q_tminus1_0treated = zeros(size(estimate_dtminus1_treated));
for z=1:size(estimate_dtminus1_treated)
    %fprintf('Within iterations No %i',z);
 estimate_untreateddiff_treated1=@(delta) estimate_untreateddiff_treated_cdf_onlypro(gammahat,x1,x2,x3,x4,D,DY2,prob,delta)-estimate_dtminus1_treated(p,1);
q_tminus1_0treated(z,1) = fzero(estimate_untreateddiff_treated1,1);
end
estimate_untreated_treated_q = @(y_input) estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input) - tau;
q_0 = fzero(estimate_untreated_treated_q,1);


iqtt(k,1) = q_1-q_0;
end
qtt_mean = mean(iqtt,1);

qtt_std = std(iqtt);

%Callaway and Li no normalization

N=1000;
T=3;
nloop = 200;
rng(1);
tau = 0.5;
x0_full = ones(N,nloop);
x1_full = normrnd(0,1,[N,nloop]);
x2_full = normrnd(0,2,[N,nloop]);
x3_full = normrnd(0,3,[N,nloop]);
x4_full = normrnd(0,4,[N,nloop]);

iqtt = zeros(nloop,1);
%both first stage estimates are correct
for k=1:nloop
fprintf('Iterations No %i ',k);
x0 = x0_full(:,k);
x1 =x1_full(:,k);
x2 =x2_full(:,k);
x3 =x3_full(:,k);
x4 =x4_full(:,k);

pro = exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4)./(1+exp(-0.25*x1 - 0.5*x2 - 0.75*x3 -1*x4));

D=zeros(size(x1));
for p=1:size(D)
D(p,1)=binornd(1,pro(p,1));
end
eta = zeros(size(D));
for p = 1:size(D)
    if D(p,1) ==0
eta(p,1)= normrnd(0,1); 
    else
    eta(p,1)= normrnd(1,1);
    end
end
v = normrnd(0,1,[N,T]);
Y1 = 0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,1);
Y2 =1 +0.50*x1 + 0.75*x2 +1.00*x3 + 1.5*x4 + eta + v(:,2);
Y3_1 =2*(0.75*x1 + 0.5*x2 +0.75*x3 + 1*x4) + eta + v(:,3);
Y3_0 =(2 +0.25*x1 + 0.5*x2 +0.75*x3 + 1*x4 + eta + v(:,3));
Y3= D.*Y3_1 +(1-D).*Y3_0;
DY1 = Y2 -Y1;
DY2 = Y3-Y2;

Y3_0_1 = zeros(size(D));
for p=1:size(Y3_0_1)
    if D(p,1)==0
        Y3_0_1 (p,1) = 0;
    else
       Y3_0_1 (p,1) = Y3_0(p,1); 
    end
end
Y3_01 = Y3_0_1(Y3_0_1~=0);

DY1_1 = zeros(size(D));
for p=1:size(DY1_1)
    if D(p,1) ==0
        DY1_1(p,1) =0;
    else
        DY1_1(p,1) = DY1(p,1);
    end
end
DY11 = DY1_1(DY1_1~=0);

DY2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_0(p,1) =DY2(p,1);
    else
        DY2_0(p,1) = 0;
    end
end
DY20 = DY2_0(DY2_0~=0);

DY2_1 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        DY2_1(p,1) =0;
    else
        DY2_1(p,1) =DY2(p,1);
    end
end
DY21 = DY2_1(DY2_1~=0);

x0_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x0_0(p,1) =x0(p,1);
    else
        x0_0(p,1) = 0;
    end
end
x00 = x0_0(x0_0~=0);

x1_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x1_0(p,1) =x1(p,1);
    else
        x1_0(p,1) = 0;
    end
end
x10 = x1_0(x1_0~=0);

x2_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x2_0(p,1) =x2(p,1);
    else
        x2_0(p,1) = 0;
    end
end
x20 = x2_0(x2_0~=0);

x3_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x3_0(p,1) =x3(p,1);
    else
        x3_0(p,1) = 0;
    end
end
x30 = x3_0(x3_0~=0);

x4_0 =zeros(size(D));
for p=1:size(D)
    if D(p,1)==0
        x4_0(p,1) =x4(p,1);
    else
        x4_0(p,1) = 0;
    end
end
x40 = x4_0(x4_0~=0);

Y1_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y1_treat(p,1) =0;
    else
        Y1_treat(p,1) =Y1(p,1);
    end
end
Y1_treat1=Y1_treat(Y1_treat~=0);

Y2_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y2_treat(p,1) =0;
    else
        Y2_treat(p,1) =Y2(p,1);
    end
end
Y2_treat1=Y2_treat(Y2_treat~=0);

Y3_treat = zeros(size(D));
for p =1:size(D)
    if D(p,1) ==0
        Y3_treat(p,1) =0;
    else
        Y3_treat(p,1) =Y3(p,1);
    end
end
Y3_treat1=Y3_treat(Y3_treat~=0);





%First Stage: Estimate beta for conditional cdf


%First Stage: Estimate logit
%initial = [0,-0.25,-0.5,-0.75,-1]; %Estimate gammas
%W = eye(5);
%options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
%[gammahat] = fminsearch(@(gamma) GMM_pro(gamma,x0,x1,x2,x3,x4,D,W),initial,options);
mynegloglik = @(gamma) -sum(log(binopdf(D,1,exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)./(1+exp(gamma(1)*x1 + gamma(2)*x2 + gamma(3)*x3 +gamma(4)*x4)))));
initial = [-0.25,-0.5,-0.75,-1];
opts = optimset('fminsearch');
opts.MaxFunEvals = Inf;
opts.MaxIter = 10000;
gammahat = fminsearch(mynegloglik,initial,opts);


%Generate empirical cdfs and cdf estimates
[f,x_y3treat] = ecdf(Y3_treat1);
q_1 = interp1(f,x_y3treat,tau,'next');

prob =mean(D);
%delta=0 This is a value to test the function below. Also, generating
%doubly robust portion
%estimate_untreateddiff_treated=estimate_untreateddiff_treated_cdf(gammahat,betahat,x0,x1,x2,x3,x4,D,DY2,prob,delta);



%estimate_tminus2_treated = estimate_treated_tminus2_cdf(Y1_treat1,delta1);

%estimate_dtminus1_treated = estimate_treated_dtminus1_cdf(DY21,delta2);

estimate_tminus2_treated =zeros(size(Y1_treat1));
%generating empirical cdf at values of Y1_treat1 (using the first period data of the treated) 
for p=1:size(Y1_treat1)
    indicator_Y1treat =double((Y1_treat1 <= Y1_treat1(p,1)));
    estimate_tminus2_treated(p,1)= mean(indicator_Y1treat,1);
end

[f1,x_y2treat] = ecdf(Y2_treat1);
q_tminus1 = zeros(size(Y2_treat1));
for p=1:size(Y2_treat1)
    q_tminus1(p,1) = interp1(f1,x_y2treat,estimate_tminus2_treated(p,1),'next');
end

estimate_dtminus1_treated = zeros(size(DY11));
for p=1:size(DY11)
    indicator_DY11 = double((DY11 <= DY11(p,1)));
    estimate_dtminus1_treated(p,1)=mean(indicator_DY11,1);
end

x_deltaY0t = -12:.001:12;
x_deltaY0t =x_deltaY0t';

%estimate_untreateddiff_treated=zeros(size(estimate_dtminus1_treated));
q_tminus1_0treated = zeros(size(estimate_dtminus1_treated));
for z=1:size(estimate_dtminus1_treated)
    %fprintf('Within iterations No %i',z);
 estimate_untreateddiff_treated1=@(delta) estimate_untreateddiff_treated_cdf_nonorm(gammahat,x1,x2,x3,x4,D,DY2,prob,delta)-estimate_dtminus1_treated(p,1);
q_tminus1_0treated(z,1) = fzero(estimate_untreateddiff_treated1,1);
end
estimate_untreated_treated_q = @(y_input) estimate_untreated_treated_cdf(q_tminus1_0treated,q_tminus1,y_input) - tau;
q_0 = fzero(estimate_untreated_treated_q,1);


iqtt(k,1) = q_1-q_0;
end
qtt_mean = mean(iqtt,1);

qtt_std = std(iqtt);



%Test data to compare to estimate at the end
N_temp=1000000;
T_temp=3;
rng(1);
x0_temp = ones(N_temp,1);
x1_temp = normrnd(0,1,[N_temp,1]);
x2_temp = normrnd(0,2,[N_temp,1]);
x3_temp = normrnd(0,3,[N_temp,1]);
x4_temp = normrnd(0,4,[N_temp,1]);

pro_temp = exp(-0.25*x1_temp - 0.5*x2_temp - 0.75*x3_temp -1*x4_temp)./(1+exp(-0.25*x1_temp - 0.5*x2_temp - 0.75*x3_temp -1*x4_temp));

D_temp=zeros(size(x1_temp));
for p=1:size(D_temp)
D_temp(p,1)=binornd(1,pro_temp(p,1));
end
eta_temp = zeros(size(D_temp));
for p = 1:size(D_temp)
    if D_temp(p,1) ==0
eta_temp(p,1)= normrnd(0,1); 
    else
    eta_temp(p,1)= normrnd(1,1);
    end
end
v_temp = normrnd(0,1,[N_temp,T_temp]);
Y1_temp = 0.25*x1_temp + 0.5*x2_temp +0.75*x3_temp + 1*x4_temp + eta_temp + v_temp(:,1);
Y2_temp =1 +0.50*x1_temp + 0.75*x2_temp +1.00*x3_temp + 1.5*x4_temp + eta_temp + v_temp(:,2);
Y3_1_temp =2*(0.75*x1_temp + 0.5*x2_temp +0.75*x3_temp + 1*x4_temp) + eta_temp + v_temp(:,3);
Y3_0_temp =(2 +0.25*x1_temp + 0.5*x2_temp +0.75*x3_temp + 1*x4_temp + eta_temp + v_temp(:,3));
Y3_temp= D_temp.*Y3_1_temp +(1-D_temp).*Y3_0_temp;
DY1_temp = Y2_temp -Y1_temp;
DY2_temp = Y3_temp-Y2_temp;


Y3_0_1_temp = zeros(size(D_temp));
for p=1:size(Y3_0_1_temp)
    if D_temp(p,1)==0
        Y3_0_1_temp (p,1) = 0;
    else
       Y3_0_1_temp (p,1) = Y3_0_temp(p,1); 
    end
end
Y3_01_temp = Y3_0_1_temp(Y3_0_1_temp~=0);


Y3_treat_temp = zeros(size(D_temp));
for p =1:size(D_temp)
    if D_temp(p,1) ==0
        Y3_treat_temp(p,1) =0;
    else
        Y3_treat_temp(p,1) =Y3_temp(p,1);
    end
end
Y3_treat1_temp=Y3_treat_temp(Y3_treat_temp~=0);


Q_1 = quantile(Y3_treat1_temp,0.75);

Q_0 = quantile(Y3_01_temp,0.75);

truth = Q_1-Q_0;
average_bias = qtt_mean -truth;
rmse_component = iqtt-truth;


RMSE = sqrt(sum((rmse_component .^2),1)/nloop);

iqtt = rmmissing(iqtt);