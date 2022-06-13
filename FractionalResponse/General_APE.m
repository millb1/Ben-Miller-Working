N=1000;
T=4;
nloop = 500;
rng(1);
pie =3.14159;

x_full = normrnd(0,4,[N,T,nloop]);
x25 = prctile(x_full(:,1,:),25,1);
x50 = prctile(x_full(:,2,:),50,1);
x75 = prctile(x_full(:,3,:),75,1);
x90 = prctile(x_full(:,4,:),90,1);


z_full = rand([N,T,nloop]);
z25 = prctile(z_full(:,1,:),25,1);
z50 = prctile(z_full(:,2,:),50,1);
z75 = prctile(z_full(:,3,:),75,1);
z90 = prctile(z_full(:,4,:),90,1);

ut_full = normrnd(0,1,[N,T,nloop]);
v_full = normrnd(0,1,[N,T,nloop]);
%v_full = random('Logistic',0,1,[N,T,nloop]);
%v_full = chi2rnd(1,[N,T,nloop]);
%delz = 1;
%delx = 0;
%Td = z*delz + x*delx + ut;         %Thing inside the t
%t = double((Td >= 0));             %Making the indicator for the thing inside t

xbar = mean(x_full,2);                  %Mean of x of firm over T years
xbar25 = prctile(xbar,25);
xbar50 = prctile(xbar,50);
xbar75 = prctile(xbar,75);
xbar90 = prctile(xbar,90);

zbar = mean(z_full,2);
zbar25 = prctile(zbar,25);
zbar50 = prctile(zbar,50);
zbar75 = prctile(zbar,75);
zbar90 = prctile(zbar,90);

%tbar = ones([1,100]);

x0star = squeeze([x25(1,1,:),x50(1,1,:),x75(1,1,:),x90(1,1,:)])';
z0star = squeeze([z25(1,1,:),z50(1,1,:),z75(1,1,:),z90(1,1,:)])';

%getting partial effects for each iternation
idelLdelx1 = ones([nloop,4]);
idelLdelx2 = ones([nloop,4]);
idelLdelh11 = ones([nloop,4]);
idelLdelh21 = ones([nloop,4]);
idelLdelh31 = ones([nloop,4]);
idelLdelh12 = ones([nloop,4]);
idelLdelh22 = ones([nloop,4]);
idelLdelh32 = ones([nloop,4]);
iSigmal1_mean = ones([nloop,4]);
iSigmal2_mean = ones([nloop,4]);
idelLdelx1c = ones([nloop,4]);
idelLdelx2c = ones([nloop,4]);
iSigmal1c_mean = ones([nloop,4]);
iSigmal2c_mean = ones([nloop,4]);

idelLdelx1a = ones([nloop,4]);
idelLdelx2a = ones([nloop,4]);
idelLdelh11a = ones([nloop,4]);
idelLdelh21a = ones([nloop,4]);
idelLdelh31a = ones([nloop,4]);
idelLdelh12a = ones([nloop,4]);
idelLdelh22a = ones([nloop,4]);
idelLdelh32a = ones([nloop,4]);
iSigmal1a_mean = ones([nloop,4]);
iSigmal2a_mean = ones([nloop,4]);
idelLdelx1ac = ones([nloop,4]);
idelLdelx2ac = ones([nloop,4]);
iSigmal1ac_mean = ones([nloop,4]);
iSigmal2ac_mean = ones([nloop,4]);

%Parameters
delz = 1;
delx = 0;
pi = [1,1,1]';
beta1 = 1;
beta2 = 2;
alpha1 = 1;
alpha2 = 2;
lambda1 = 1;
lambda2 = 2;
zeta1 = 1;
zeta2 = 1;

for k = 1:nloop
fprintf('Iterations No %i ',k);
x = x_full(:,:,k);
z = z_full(:,:,k);
ut = ut_full(:,:,k);
v = v_full(:,:,k);
    
Td = z*delz + x*delx + ut;    
t = double((Td >= 0));
    
h = [xbar(:,:,k),zbar(:,:,k),mean(t,2)];   

r1 = zeta1*ut + v;
r2 = zeta2*ut + v;

G1 = (exp(x.*beta1+t*alpha1+lambda1*h*pi+r1))./(1 + exp(x.*beta1+t.*alpha1+lambda1*h*pi+r1) + exp(x*beta2+t*alpha2+lambda2*h*pi+r2));
G2 = (exp(x*beta2+t*alpha2+lambda2*h*pi+r2))./(1 + exp(x*beta1+t*alpha1+lambda1*h*pi+r1) + exp(x*beta2+t*alpha2+lambda2*h*pi+r2));
G3 = 1 - G1 - G2;
G1(G3 <0 ) =round(G1(G3<0));
G2(G3<0) =round(G2(G3<0));
G3 = 1 - G1 - G2;

outcomes = zeros(N,T,3);
for i = 1:N
    for j = 1:T
        p = [G1(i,j),G2(i,j),G3(i,j)];
        y = mnrnd(100,p,1);
        outcomes(i,j,1) = y(1);
        outcomes(i,j,2) = y(2);
        outcomes(i,j,3) = y(3);
    end
end
            
fracres = outcomes./100;
      
%% GMM
%Gauss-Legauerre quadrature setup
node = 10;
syms xx
c = double(coeffs(laguerreL(node,xx),'All'));
xx= sort(roots(c));
o = ones([1,node]);
for i = 1:node
    o(i) = xx(i)/(((node+1)*laguerreL(node+1,xx(i)))^2);
end
u = ones([N,T,node]);
for i = 1:node
    u(:,:,i) = xx(i).*u(:,:,i);
end
w = ones([N,T,node]);
for i = 1:node
    w(:,:,i) = o(i).* w(:,:,i);
end 
%First Stage: Estimate Delta
initial = zeros(2,1)'; %[delz,delx]
W = eye(2);
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter', 10000 );
[deltahat] = fminsearch(@(delta) GMM_delta_reduced(delta,z,x,t,W),initial,options);

%Second Stage: Estimate Theta
initial = [1,2,1,2,1,2,1,1,1,1,1]; %[beta1,beta2,alpha1,alpha2,lambda1,lambda2,pi(1),pi(2),pi(4),zeta1,zeta2];
W= eye(11);
[thetahat] = fminsearch(@(theta) GMM_theta_reduced_glquad(theta,deltahat,z,x,h,u,t,w,fracres,W,pie),initial,options);    
    
%Estimated Parameters
delzhat = deltahat(1);
delxhat = deltahat(2);
beta1hat = thetahat(1);
beta2hat = thetahat(2);
alpha1hat = thetahat(3);
alpha2hat = thetahat(4);
lambda1hat = thetahat(5);
lambda2hat = thetahat(6);
pi1hat = thetahat(7);
pi2hat = thetahat(8);
pi3hat = thetahat(9);
zeta1hat = thetahat(10);
zeta2hat = thetahat(11);

%% estimated APE    
    
x0 = x0star(k,:);
z0 = z0star(k,:);
t0 = ones([1,4]);
h0 = [xbar25(k),xbar50(k),xbar75(k),xbar90(k); zbar25(k),zbar50(k),zbar75(k),zbar90(k); t0]';

dit1_0 = beta1hat.*x0 + alpha1hat.*t0 + lambda1hat.*pi1hat.*h0(:,1)' + lambda1hat.*pi2hat.*h0(:,2)' + lambda1hat.*pi3hat.*h0(:,3)';
dit2_0 = beta2hat.*x0 + alpha2hat.*t0 + lambda2hat.*pi1hat.*h0(:,1)' + lambda2hat.*pi2hat.*h0(:,2)' + lambda2hat.*pi3hat.*h0(:,3)';

rng(1)
sims = 1000;
u2 = normrnd(0,1,[sims,1]);

Lprime1_0 = (exp(dit1_0 + zeta1hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit1_0+zeta1hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;
Lprime2_0 = (exp(dit2_0 + zeta2hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit2_0+zeta2hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;

delLdelx1 = mean(beta1hat*Lprime1_0,1);
idelLdelx1(k,:) = delLdelx1;

%delLdelx1_std = sqrt(mean(((beta1hat*Lprime1_0)-delLdelx1).^2));
delLdelx2 = mean(beta2hat*Lprime2_0,1);
%delLdelx2_std=sqrt(mean(((beta2hat*Lprime2_0)-delLdelx2).^2));
idelLdelx2(k,:) = delLdelx2;

delLdelh11 = mean(lambda1hat*pi1hat*Lprime1_0,1);
%delLdelh11_std = sqrt(mean(((lambda1hat*pi1hat*Lprime1_0)-delLdelh11).^2));
idelLdelh11(k,:) = delLdelh11;

delLdelh21 = mean(lambda1hat*pi2hat*Lprime1_0,1);
%delLdelh21_std = sqrt(mean(((lambda1hat*pi2hat*Lprime1_0)-delLdelh21).^2));
idelLdelh21(k,:) = delLdelh21;

delLdelh31 = mean(lambda1hat*pi3hat*Lprime1_0,1);
%delLdelh31_std = sqrt(mean(((lambda1hat*pi3hat*Lprime1_0)-delLdelh31).^2));
idelLdelh31(k,:) = delLdelh31;

delLdelh12 = mean(lambda2hat*pi1hat*Lprime2_0,1);
%delLdelh12_std = sqrt(mean(((lambda2hat*pi1hat*Lprime2_0)-delLdelh12).^2));
idelLdelh12(k,:) = delLdelh12;

delLdelh22 = mean(lambda2hat*pi2hat*Lprime2_0,1);
%delLdelh22_std = sqrt(mean(((lambda2hat*pi2hat*Lprime2_0)-delLdelh22).^2));
idelLdelh22(k,:) = delLdelh22;

delLdelh32 = mean(lambda2hat*pi3hat*Lprime2_0,1);
%delLdelh32_std = sqrt(mean(((lambda2hat*pi3hat*Lprime2_0)-delLdelh32).^2));
idelLdelh32(k,:) = delLdelh32;        

% APE ts
ditt11 = beta1hat.*x0  + alpha1hat + lambda1hat.*pi1hat.*h0(:,1)' + lambda1hat.*pi2hat.*h0(:,2)' + lambda1hat.*pi3hat.*h0(:,3)';
ditt12 = beta2hat.*x0  + alpha2hat + lambda2hat.*pi1hat.*h0(:,1)' + lambda2hat.*pi2hat.*h0(:,2)' + lambda2hat.*pi3hat.*h0(:,3)';

ditt01 = beta1hat.*x0  + lambda1hat.*pi1hat.*h0(:,1)' + lambda1hat.*pi2hat.*h0(:,2)' + lambda1hat.*pi3hat.*h0(:,3)';
ditt02 = beta2hat.*x0  + lambda2hat.*pi1hat.*h0(:,1)' + lambda2hat.*pi2hat.*h0(:,2)' + lambda2hat.*pi3hat.*h0(:,3)';

L1t0 = exp(ditt01+zeta1hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L1t1 = exp(ditt11+zeta1hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

L2t0 = exp(ditt02+zeta2hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L2t1 = exp(ditt12+zeta2hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

Sigmal1_mean = mean(mean(L1t1-L1t0,1),1);
%Sigmal1_std=sqrt(mean(((L1t1-L1t0)-Sigmal1_mean).^2));
iSigmal1_mean(k,:) = Sigmal1_mean;
Sigmal2_mean = mean(mean(L2t1-L2t0,1),1);
%Sigmal2_std=sqrt(mean(((L2t1-L2t0)-Sigmal2_mean).^2));
iSigmal2_mean(k,:) = Sigmal2_mean;

%APE true
dit1_a = beta1.*x0  + alpha1.*t0 + lambda1.*1.*h0(:,1)' + lambda1.*1.*h0(:,2)' + lambda1.*1.*h0(:,3)';
dit2_a = beta2.*x0  + alpha2.*t0 + lambda2.*1.*h0(:,1)' + lambda2.*1.*h0(:,2)' + lambda2.*1.*h0(:,3)';   
%"actual" continuous APEs 
%v2 = random('Logistic',0,1,[sims,1]);
v2 = normrnd(0,1,[sims,1]);
%v2 = chi2rnd(1,[sims,1]);

Lprime1_a = (exp(dit1_a + zeta1*u2 + v2).*(1 + exp(dit1_a+zeta1*u2 + v2) + exp(dit2_a + zeta2*u2 + v2)) - exp(2*(dit1_a+zeta1*u2 + v2)))./(1 + exp(dit1_a+zeta1*u2+v2) + exp(dit2_a + zeta2*u2+v2)).^2;
Lprime2_a = (exp(dit2_a + zeta2*u2 + v2).*(1 + exp(dit1_a+zeta1*u2 + v2) + exp(dit2_a + zeta2*u2 + v2)) - exp(2*(dit2_a+zeta2*u2 + v2)))./(1 + exp(dit1_a+zeta1*u2+v2) + exp(dit2_a + zeta2*u2+v2)).^2;

delLdelx1a = mean(beta1*Lprime1_a,1);
idelLdelx1a(k,:) = delLdelx1a;

delLdelx2a = mean(beta2*Lprime2_a,1);
idelLdelx2a(k,:) = delLdelx2a;
  
delLdelh11a = mean(lambda1*1*Lprime1_a,1);
idelLdelh11a(k,:) = delLdelh11a;

delLdelh21a = mean(lambda1*1*Lprime1_a,1);
idelLdelh21a(k,:) = delLdelh21a;

delLdelh31a = mean(lambda1*1*Lprime1_a,1);
idelLdelh31a(k,:) = delLdelh31a;

delLdelh12a = mean(lambda2*1*Lprime2_a,1);
idelLdelh12a(k,:) = delLdelh12a;

delLdelh22a = mean(lambda2*1*Lprime2_a,1);
idelLdelh22a(k,:) = delLdelh22a;

delLdelh32a = mean(lambda2*1*Lprime2_a,1);
idelLdelh32a(k,:) = delLdelh32a;

%"actual ds with varying values of t
ditt11_a = beta1.*x0  + alpha1 + lambda1.*1.*h0(:,1)' + lambda1.*1.*h0(:,2)' + lambda1.*1.*h0(:,3)';
ditt12_a = beta2.*x0  + alpha2 + lambda2.*1.*h0(:,1)' + lambda2.*1.*h0(:,2)' + lambda2.*1.*h0(:,3)';

ditt01_a = beta1.*x0  + lambda1.*1.*h0(:,1)' + lambda1.*1.*h0(:,2)' + lambda1.*1.*h0(:,3)';
ditt02_a = beta2.*x0  + lambda2.*1.*h0(:,1)' + lambda2.*1.*h0(:,2)' + lambda2.*1.*h0(:,3)';


% "actual" discrete APEs 
L1t0a = exp(ditt01_a+zeta1*u2+v2)./(1 + exp(ditt01_a+zeta1*u2+v2) + exp(ditt02_a + zeta2*u2+v2));
L1t1a = exp(ditt11_a+zeta1*u2+v2)./(1 + exp(ditt11_a+zeta1*u2+v2) + exp(ditt12_a + zeta2*u2+v2));

L2t0a = exp(ditt02_a+zeta2*u2+v2)./(1 + exp(ditt01_a+zeta1*u2+v2) + exp(ditt02_a + zeta2*u2+v2));
L2t1a = exp(ditt12_a+zeta2*u2+v2)./(1 + exp(ditt11_a+zeta1*u2+v2) + exp(ditt12_a + zeta2*u2+v2));

Sigmal1a_mean =  mean(mean(L1t1a-L1t0a,1),1);
iSigmal1a_mean(k,:) = Sigmal1a_mean;

Sigmal2a_mean = mean(mean(L2t1a-L2t0a,1),1);
iSigmal2a_mean(k,:) = Sigmal2a_mean;

%APE averaged over c and e
%Estimated APE

ndraws = 300;
idelLdelx1c_draws = ones([ndraws,4]);
idelLdelx2c_draws = ones([ndraws,4]);
iSigmal1c_mean_draws = ones([ndraws,4]);
iSigmal2c_mean_draws = ones([ndraws,4]);
idelLdelx1ac_draws = ones([ndraws,4]);
idelLdelx2ac_draws = ones([ndraws,4]);
iSigmal1ac_mean_draws = ones([ndraws,4]);
iSigmal2ac_mean_draws = ones([ndraws,4]);


for j=1:N
dit1_c = beta1hat.*x0 + alpha1hat.*t0 + lambda1hat.*pi1hat.*h(j,1)' + lambda1hat.*pi2hat.*h(j,2)' + lambda1hat.*pi3hat.*h(j,3)';
dit2_c = beta2hat.*x0 + alpha2hat.*t0 + lambda2hat.*pi1hat.*h(j,1)' + lambda2hat.*pi2hat.*h(j,2)' + lambda2hat.*pi3hat.*h(j,3)';

sims = 1000;
u3 = normrnd(0,1,[ndraws,1]);

Lprime1_c = (exp(dit1_c + zeta1hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit1_c+zeta1hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
Lprime2_c = (exp(dit2_c + zeta2hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit2_c+zeta2hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
delLdelx1c = mean(beta1hat*Lprime1_c,1);
idelLdelx1c_draws(j,:) = delLdelx1c;
delLdelx2c = mean(beta2hat*Lprime2_c,1);
idelLdelx2c_draws(j,:) = delLdelx2c;

ditt11c = beta1hat.*x0  + alpha1hat + lambda1hat.*pi1hat.*h(j,1)' + lambda1hat.*pi2hat.*h(j,2)' + lambda1hat.*pi3hat.*h(j,3)';
ditt12c = beta2hat.*x0  + alpha2hat + lambda2hat.*pi1hat.*h(j,1)' + lambda2hat.*pi2hat.*h(j,2)' + lambda2hat.*pi3hat.*h(j,3)';

ditt01c = beta1hat.*x0  + lambda1hat.*pi1hat.*h(j,1)' + lambda1hat.*pi2hat.*h(j,2)' + lambda1hat.*pi3hat.*h(j,3)';
ditt02c = beta2hat.*x0  + lambda2hat.*pi1hat.*h(j,1)' + lambda2hat.*pi2hat.*h(j,2)' + lambda2hat.*pi3hat.*h(j,3)';

L1t0c = exp(ditt01c+zeta1hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L1t1c = exp(ditt11c+zeta1hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

L2t0c = exp(ditt02c+zeta2hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L2t1c = exp(ditt12c+zeta2hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

Sigmal1_meanc = mean(mean(L1t1c-L1t0c,1),1);
iSigmal1c_mean_draws(j,:) = Sigmal1_meanc;
Sigmal2_meanc = mean(mean(L2t1c-L2t0c,1),1);
iSigmal2c_mean_draws(j,:) = Sigmal2_meanc;

%Actual APE

idelLdelx1ac_draws_v = ones([ndraws,4]);
idelLdelx2ac_draws_v = ones([ndraws,4]);
iSigmal1ac_mean_draws_v = ones([ndraws,4]);
iSigmal2ac_mean_draws_v = ones([ndraws,4]);



for g=1:ndraws
   
 v3 = normrnd(0,1,[ndraws,1]); 
dit1_ca = beta1.*x0 + alpha1.*t0 + lambda1.*1.*h(j,1)' + lambda1.*1.*h(j,2)' + lambda1.*1.*h(j,3)';
dit2_ca = beta2.*x0 + alpha2.*t0 + lambda2.*1.*h(j,1)' + lambda2.*1.*h(j,2)' + lambda2.*1.*h(j,3)';

Lprime1_ca = (exp(dit1_ca + zeta1*u3(g,:) + v3).*(1 + exp(dit1_ca+zeta1*u3(g,:) + v3) + exp(dit2_ca + zeta2*u3(g,:) + v3)) - exp(2*(dit1_ca+zeta1*u3(g,:) + v3)))./(1 + exp(dit1_ca+zeta1*u3(g,:)+v3) + exp(dit2_ca + zeta2*u3(g,:)+v3)).^2;
Lprime2_ca = (exp(dit2_ca + zeta2*u3(g,:) + v3).*(1 + exp(dit1_ca+zeta1*u3(g,:) + v3) + exp(dit2_ca + zeta2*u3(g,:) + v3)) - exp(2*(dit2_ca+zeta2*u3(g,:) + v3)))./(1 + exp(dit1_ca+zeta1*u3(g,:)+v3) + exp(dit2_ca + zeta2*u3(g,:)+v3)).^2;
 
ditt11ca = beta1hat.*x0  + alpha1 + lambda1.*1.*h(j,1)' + lambda1.*1.*h(j,2)' + lambda1.*1.*h(j,3)';
ditt12ca = beta2hat.*x0  + alpha2 + lambda2.*1.*h(j,1)' + lambda2.*1.*h(j,2)' + lambda2.*1.*h(j,3)';

ditt01ca = beta1.*x0  + lambda1.*1.*h(j,1)' + lambda1.*1.*h(j,2)' + lambda1.*1.*h(j,3)';
ditt02ca = beta2.*x0  + lambda2.*1.*h(j,1)' + lambda2.*1.*h(j,2)' + lambda2.*1.*h(j,3)';

L1t0ca = exp(ditt01ca+zeta1*u3(g,:)+v3)./(1 + exp(ditt01ca+zeta1*u3(g,:)+v3) + exp(ditt02ca + zeta2*u3(g,:)+v3));
L1t1ca = exp(ditt11ca+zeta1*u3(g,:)+v3)./(1 + exp(ditt11ca+zeta1*u3(g,:)+v3) + exp(ditt12ca + zeta2*u3(g,:)+v3));

L2t0ca = exp(ditt02ca+zeta2*u3(g,:)+v3)./(1 + exp(ditt01ca+zeta1*u3(g,:)+v3) + exp(ditt02ca + zeta2*u3(g,:)+v3));
L2t1ca = exp(ditt12ca+zeta2*u3(g,:)+v3)./(1 + exp(ditt11ca+zeta1*u3(g,:)+v3) + exp(ditt12ca + zeta2*u3(g,:)+v3));

delLdelx1ca =mean(beta1*Lprime1_ca,1);
delLdelx2ca =mean(beta2*Lprime2_ca,1);
Sigmal1_meanca = mean(mean(L1t1ca-L1t0ca,1),1);
Sigmal2_meanca = mean(mean(L2t1ca-L2t0ca,1),1);

idelLdelx1ac_draws_v(g,:) = delLdelx1ca;
idelLdelx2ac_draws_v(g,:) = delLdelx2ca;
iSigmal1ac_mean_draws_v(g,:) = Sigmal1_meanca;
iSigmal2ac_mean_draws_v(g,:) = Sigmal2_meanca;

end
idelLdelx1ac_draws(j,:) = mean(idelLdelx1ac_draws_v,1);
idelLdelx2ac_draws(j,:) = mean(idelLdelx2ac_draws_v,1);
iSigmal1ac_mean_draws(j,:) = mean(iSigmal1ac_mean_draws_v,1);
iSigmal2ac_mean_draws(j,:) = mean(iSigmal2ac_mean_draws_v,1);

end
idelLdelx1c(k,:) = mean(idelLdelx1c_draws,1);
idelLdelx2c(k,:) = mean(idelLdelx2c_draws,1);
iSigmal1c_mean(k,:) = mean(iSigmal1c_mean_draws,1);
iSigmal2c_mean(k,:) = mean(iSigmal2c_mean_draws,1);

idelLdelx1ac(k,:) = mean(idelLdelx1ac_draws,1);
idelLdelx2ac(k,:) = mean(idelLdelx2ac_draws,1);
iSigmal1ac_mean(k,:) = mean(iSigmal1ac_mean_draws,1);
iSigmal2ac_mean(k,:) = mean(iSigmal2ac_mean_draws,1);



end


%Averages
meaneAPEx1 = mean(idelLdelx1); 
sdeAPEx1 = sqrt(mean((idelLdelx1 - meaneAPEx1).^2));

meaneAPEx2 = mean(idelLdelx2); 
sdeAPEx2 = sqrt(mean((idelLdelx2 - meaneAPEx2).^2));


meaneAPEh11 = mean(idelLdelh11);
sdeAPEh11 = sqrt(mean((idelLdelh11 - meaneAPEh11).^2));


meaneAPEh21 = mean(idelLdelh21);
sdeAPEh21 = sqrt(mean((idelLdelh21 - meaneAPEh21).^2));


meaneAPEh31 = mean(idelLdelh31);
sdeAPEh31 = sqrt(mean((idelLdelh31 - meaneAPEh31).^2));


meaneAPEh12 = mean(idelLdelh12);
sdeAPEh12 = sqrt(mean((idelLdelh12 - meaneAPEh12).^2));


meaneAPEh22 = mean(idelLdelh22); 
sdeAPEh22 = sqrt(mean((idelLdelh22 - meaneAPEh22).^2));


meaneAPEh32 = mean(idelLdelh32);
sdeAPEh32 = sqrt(mean((idelLdelh32 - meaneAPEh32).^2));


meaneAPEt1 = mean(iSigmal1_mean); 
sdeAPEt1 = sqrt(mean((iSigmal1_mean - meaneAPEt1).^2));

meaneAPEt2 = mean(iSigmal2_mean); 
sdeAPEt2 = sqrt(mean((iSigmal2_mean - meaneAPEt2).^2));

meaneAPEx1c = mean(idelLdelx1c);
sdeAPEx1c = sqrt(mean((idelLdelx1c - meaneAPEx1c).^2));

meaneAPEx2c = mean(idelLdelx2c);
sdeAPEx2c = sqrt(mean((idelLdelx2c - meaneAPEx2c).^2));

meaneAPEt1c = mean(iSigmal1c_mean); 
sdeAPEt1c = sqrt(mean((iSigmal1c_mean - meaneAPEt1c).^2));

meaneAPEt2c = mean(iSigmal2c_mean); 
sdeAPEt2c = sqrt(mean((iSigmal2c_mean - meaneAPEt2c).^2));

%Truth
meaneAPEx1a = mean(idelLdelx1a);

meaneAPEx2a = mean(idelLdelx2a);

meaneAPEh11a = mean(idelLdelh11a);

meaneAPEh21a = mean(idelLdelh21a);

meaneAPEh31a = mean(idelLdelh31a);

meaneAPEh12a = mean(idelLdelh12a);

meaneAPEh22a = mean(idelLdelh22a);

meaneAPEh32a = mean(idelLdelh32a);

meaneAPEt1a = mean(iSigmal1a_mean);

meaneAPEt2a = mean(iSigmal2a_mean);

meaneAPEx1ac = mean(idelLdelx1ac);

meaneAPEx2ac = mean(idelLdelx2ac);

meaneAPEt1ac = mean(iSigmal1ac_mean);

meaneAPEt2ac = mean(iSigmal2ac_mean);
