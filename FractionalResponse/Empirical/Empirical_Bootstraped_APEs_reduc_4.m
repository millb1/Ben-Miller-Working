


%Bootstrap Standard Errors with Intercept adjusting for endogeneity
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter
    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    z_samp = z(s_index(:,:,i),:);
    z_samp = reshape(z_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    u_temp = u(s_index(:,:,i),:,:);
    w_temp = w(s_index(:,:,i),:,:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    initial =[.0009818,.0085474,.0474713,-.2565403,.5223746,.0936121 ,.0106284,.0008201,-.3571223,-1.487394];
    W = eye(size(initial,1));
    [tempdeltahat] = fminsearch(@(delta) SGMM_Empiric_S1(delta,z_samp,x_samp,K(1),t_samp,W),initial,options);
    initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,1]';
    W = eye(size(initial,1));
    pie = 3.14159;
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2(theta,tempdeltahat,z_samp,x_samp,h_samp,t_samp,u_temp,w_temp,fracres_temp,W,pie,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(2),T,inter]);
b_ape_x2 = zeros([K(2),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

ndraws = 300;
idelLdelx1c_draws = ones([K(2),4,samps]);
idelLdelx2c_draws = ones([K(2),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps
 u3 = normrnd(0,1,[1,1,ndraws]);   
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c + zeta1hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit1_c+zeta1hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
Lprime2_c = (exp(dit2_c + zeta2hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit2_c+zeta2hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
delLdelx1c = beta1hat.*mean(Lprime1_c,3); 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*mean(Lprime2_c,3); 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c+zeta1hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L1t1c = exp(ditt11c+zeta1hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

L2t0c = exp(ditt02c+zeta2hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L2t1c = exp(ditt12c+zeta2hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

Sigmal1_meanc = mean(mean(L1t1c-L1t0c,1),3);
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = mean(mean(L2t1c-L2t0c,1),3);
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));


%% Bootstrap Simple APE
zp25 = prctile(z(:,1,:),25,1);
zp50 = prctile(z(:,2,:),50,1);
zp75 = prctile(z(:,3,:),75,1);
zp90 = prctile(z(:,4,:),90,1);
z0star = squeeze([zp25(1,1,:),zp50(1,1,:),zp75(1,1,:),zp90(1,1,:)])';
h0 = [x0star(1,:);x0star(2,:);x0star(3,:);x0star(6,:);x0star(7,:);z0star(:,2)';t0]';

b_sape_x1 = zeros([K(2),T,inter]);
b_sape_x2 = zeros([K(2),T,inter]);
b_sape_t1 = zeros([1,T,inter]);
b_sape_t2 = zeros([1,T,inter]);
b_sape_h = zeros([K(3),T,inter]);

for i = 1:inter
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

dit1_0 = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + (lambda1hat*h0*pihat)';
dit2_0 = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + (lambda2hat*h0*pihat)';

sims = 1000;
u2 = normrnd(0,1,[1,1,sims]);

Lprime1_0 = (exp(dit1_0 + zeta1hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit1_0+zeta1hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;
Lprime2_0 = (exp(dit2_0 + zeta2hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit2_0+zeta2hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;

%x's
b_sape_x1(:,:,i) = beta1hat.*mean(Lprime1_0,3);
b_sape_x2(:,:,i) = beta2hat.*mean(Lprime2_0,3);

%h
delLdelht = lambda1hat*reshape(pihat,1,1,[]).*mean(Lprime1_0,3);
b_sape_h(:,:,i) = reshape(mean(delLdelht,1),4,13)';

%t's
ditt11 = sum(x0star.*beta1hat,1)  + alpha1hat + (lambda1hat*h0*pihat)';
ditt12 = sum(x0star.*beta2hat,1)  + alpha2hat + (lambda2hat*h0*pihat)';

ditt01 = sum(x0star.*rbeta1hat,1)  + (lambda1hat*h0*pihat)';
ditt02 = sum(x0star.*beta2hat,1)  + (lambda2hat*h0*pihat)';

L1t0 = exp(ditt01+zeta1hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L1t1 = exp(ditt11+zeta1hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

L2t0 = exp(ditt02+zeta2hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L2t1 = exp(ditt12+zeta2hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

b_sape_t1(:,:,i) = mean(mean(L1t1-L1t0,3),1);
b_sape_t2(:,:,i) = mean(mean(L2t1-L2t0,3),1);
end

%Bootstraped Standard Errors
%x's
sape_x1_sd = sqrt(sum((b_sape_x1 - mean(b_sape_x1)).^2,3)./(inter-1));
sape_x2_sd = sqrt(sum((b_sape_x2 - mean(b_sape_x2)).^2,3)./(inter-1));

%t's
sape_t1_sd = sqrt(sum((b_sape_t1 - mean(b_sape_t1)).^2,3)./(inter-1));
sape_t2_sd = sqrt(sum((b_sape_t2 - mean(b_sape_t2)).^2,3)./(inter-1));

%h
sape_h_sd = sqrt(sum((b_sape_h - mean(b_sape_h)).^2,3)./(inter-1));

%Bootstrap with cre but no adjustment for endogeneity*****
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter
    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    
    initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
    W = eye(size(initial,1));
   
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2_no_endog(theta,x_samp,h_samp,t_samp,fracres_temp,W,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(1),T,inter]);
b_ape_x2 = zeros([K(1),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(1),i);
beta2hat = b_thetahat(K(1)+1:2*K(1),i);
pihat = b_thetahat(2*K(1)+1:2*K(1)+K(2),i);
alpha1hat = b_thetahat(length(thetahat)-3,i);
alpha2hat = b_thetahat(length(thetahat)-2,i);
lambda1hat = b_thetahat(length(thetahat)-1,i);
lambda2hat = b_thetahat(length(thetahat),i);

idelLdelx1c_draws = ones([K(1),4,samps]);
idelLdelx2c_draws = ones([K(1),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps    
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit1_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
Lprime2_c = (exp(dit2_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit2_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
delLdelx1c = beta1hat.*Lprime1_c; 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*Lprime2_c; 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c)./(1 + exp(ditt01c) + exp(ditt02c));
L1t1c = exp(ditt11c)./(1 + exp(ditt11c) + exp(ditt12c));

L2t0c = exp(ditt02c)./(1 + exp(ditt01c) + exp(ditt02c));
L2t1c = exp(ditt12c)./(1 + exp(ditt11c) + exp(ditt12c));

Sigmal1_meanc = L1t1c-L1t0c;
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = L2t1c-L2t0c;
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));




%****Male only data******************************

%Bootstrap Standard Errors with Intercept adjusting for endogeneity
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter

    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    z_samp = z(s_index(:,:,i),:);
    z_samp = reshape(z_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    u_temp = u(s_index(:,:,i),:,:);
    w_temp = w(s_index(:,:,i),:,:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    initial =[.0008149,.0073264,.0561136,-.102546,.8566147 ,.1234185,.0078028,.0080065,-.360758,-1.811018];
    W = eye(size(initial,1));
    [tempdeltahat] = fminsearch(@(delta) SGMM_Empiric_S1(delta,z_samp,x_samp,K(1),t_samp,W),initial,options);
   initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,1]';
    W = eye(size(initial,1));
    pie = 3.14159;
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2(theta,tempdeltahat,z_samp,x_samp,h_samp,t_samp,u_temp,w_temp,fracres_temp,W,pie,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(2),T,inter]);
b_ape_x2 = zeros([K(2),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

ndraws = 300;
idelLdelx1c_draws = ones([K(2),4,samps]);
idelLdelx2c_draws = ones([K(2),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps
 u3 = normrnd(0,1,[1,1,ndraws]);   
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c + zeta1hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit1_c+zeta1hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
Lprime2_c = (exp(dit2_c + zeta2hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit2_c+zeta2hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
delLdelx1c = beta1hat.*mean(Lprime1_c,3); 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*mean(Lprime2_c,3); 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c+zeta1hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L1t1c = exp(ditt11c+zeta1hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

L2t0c = exp(ditt02c+zeta2hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L2t1c = exp(ditt12c+zeta2hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

Sigmal1_meanc = mean(mean(L1t1c-L1t0c,1),3);
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = mean(mean(L2t1c-L2t0c,1),3);
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));


%% Bootstrap Simple APE
zp25 = prctile(z(:,1,:),25,1);
zp50 = prctile(z(:,2,:),50,1);
zp75 = prctile(z(:,3,:),75,1);
zp90 = prctile(z(:,4,:),90,1);
z0star = squeeze([zp25(1,1,:),zp50(1,1,:),zp75(1,1,:),zp90(1,1,:)])';
h0 = [x0star(1,:);x0star(2,:);x0star(3,:);x0star(6,:);x0star(7,:);z0star(:,2)';t0]';

b_sape_x1 = zeros([K(2),T,inter]);
b_sape_x2 = zeros([K(2),T,inter]);
b_sape_t1 = zeros([1,T,inter]);
b_sape_t2 = zeros([1,T,inter]);
b_sape_h = zeros([K(3),T,inter]);

for i = 1:inter
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

dit1_0 = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + (lambda1hat*h0*pihat)';
dit2_0 = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + (lambda2hat*h0*pihat)';

sims = 1000;
u2 = normrnd(0,1,[1,1,sims]);

Lprime1_0 = (exp(dit1_0 + zeta1hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit1_0+zeta1hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;
Lprime2_0 = (exp(dit2_0 + zeta2hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit2_0+zeta2hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;

%x's
b_sape_x1(:,:,i) = beta1hat.*mean(Lprime1_0,3);
b_sape_x2(:,:,i) = beta2hat.*mean(Lprime2_0,3);

%h
delLdelht = lambda1hat*reshape(pihat,1,1,[]).*mean(Lprime1_0,3);
b_sape_h(:,:,i) = reshape(mean(delLdelht,1),4,13)';

%t's
ditt11 = sum(x0star.*beta1hat,1)  + alpha1hat + (lambda1hat*h0*pihat)';
ditt12 = sum(x0star.*beta2hat,1)  + alpha2hat + (lambda2hat*h0*pihat)';

ditt01 = sum(x0star.*rbeta1hat,1)  + (lambda1hat*h0*pihat)';
ditt02 = sum(x0star.*beta2hat,1)  + (lambda2hat*h0*pihat)';

L1t0 = exp(ditt01+zeta1hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L1t1 = exp(ditt11+zeta1hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

L2t0 = exp(ditt02+zeta2hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L2t1 = exp(ditt12+zeta2hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

b_sape_t1(:,:,i) = mean(mean(L1t1-L1t0,3),1);
b_sape_t2(:,:,i) = mean(mean(L2t1-L2t0,3),1);
end

%Bootstraped Standard Errors
%x's
sape_x1_sd = sqrt(sum((b_sape_x1 - mean(b_sape_x1)).^2,3)./(inter-1));
sape_x2_sd = sqrt(sum((b_sape_x2 - mean(b_sape_x2)).^2,3)./(inter-1));

%t's
sape_t1_sd = sqrt(sum((b_sape_t1 - mean(b_sape_t1)).^2,3)./(inter-1));
sape_t2_sd = sqrt(sum((b_sape_t2 - mean(b_sape_t2)).^2,3)./(inter-1));

%h
sape_h_sd = sqrt(sum((b_sape_h - mean(b_sape_h)).^2,3)./(inter-1));

%Bootstrap with cre but no adjustment for endogeneity*****
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter
    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    
    initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
    W = eye(size(initial,1));
   
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2_no_endog(theta,x_samp,h_samp,t_samp,fracres_temp,W,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(1),T,inter]);
b_ape_x2 = zeros([K(1),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(1),i);
beta2hat = b_thetahat(K(1)+1:2*K(1),i);
pihat = b_thetahat(2*K(1)+1:2*K(1)+K(2),i);
alpha1hat = b_thetahat(length(thetahat)-3,i);
alpha2hat = b_thetahat(length(thetahat)-2,i);
lambda1hat = b_thetahat(length(thetahat)-1,i);
lambda2hat = b_thetahat(length(thetahat),i);

idelLdelx1c_draws = ones([K(1),4,samps]);
idelLdelx2c_draws = ones([K(1),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps    
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit1_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
Lprime2_c = (exp(dit2_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit2_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
delLdelx1c = beta1hat.*Lprime1_c; 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*Lprime2_c; 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c)./(1 + exp(ditt01c) + exp(ditt02c));
L1t1c = exp(ditt11c)./(1 + exp(ditt11c) + exp(ditt12c));

L2t0c = exp(ditt02c)./(1 + exp(ditt01c) + exp(ditt02c));
L2t1c = exp(ditt12c)./(1 + exp(ditt11c) + exp(ditt12c));

Sigmal1_meanc = L1t1c-L1t0c;
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = L2t1c-L2t0c;
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));



%****female only data******************************

%Bootstrap Standard Errors with Intercept adjusting for endogeneity
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter

    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    z_samp = z(s_index(:,:,i),:);
    z_samp = reshape(z_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    u_temp = u(s_index(:,:,i),:,:);
    w_temp = w(s_index(:,:,i),:,:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    initial =[.0011,.0069,.0418,-.3888,.3462,.1008,.0152,-.0110,-.2585,-1.118];
    W = eye(size(initial,1));
    [tempdeltahat] = fminsearch(@(delta) SGMM_Empiric_S1(delta,z_samp,x_samp,K(1),t_samp,W),initial,options);
   initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,1]';
    W = eye(size(initial,1));
    pie = 3.14159;
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2(theta,tempdeltahat,z_samp,x_samp,h_samp,t_samp,u_temp,w_temp,fracres_temp,W,pie,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(2),T,inter]);
b_ape_x2 = zeros([K(2),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

ndraws = 300;
idelLdelx1c_draws = ones([K(2),4,samps]);
idelLdelx2c_draws = ones([K(2),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps
 u3 = normrnd(0,1,[1,1,ndraws]);   
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c + zeta1hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit1_c+zeta1hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
Lprime2_c = (exp(dit2_c + zeta2hat*u3).*(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)) - exp(2*(dit2_c+zeta2hat*u3)))./(1 + exp(dit1_c+zeta1hat*u3) + exp(dit2_c + zeta2hat*u3)).^2;
delLdelx1c = beta1hat.*mean(Lprime1_c,3); 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*mean(Lprime2_c,3); 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c+zeta1hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L1t1c = exp(ditt11c+zeta1hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

L2t0c = exp(ditt02c+zeta2hat*u3)./(1 + exp(ditt01c+zeta1hat*u3) + exp(ditt02c + zeta2hat*u3));
L2t1c = exp(ditt12c+zeta2hat*u3)./(1 + exp(ditt11c+zeta1hat*u3) + exp(ditt12c + zeta2hat*u3));

Sigmal1_meanc = mean(mean(L1t1c-L1t0c,1),3);
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = mean(mean(L2t1c-L2t0c,1),3);
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));


%% Bootstrap Simple APE
zp25 = prctile(z(:,1,:),25,1);
zp50 = prctile(z(:,2,:),50,1);
zp75 = prctile(z(:,3,:),75,1);
zp90 = prctile(z(:,4,:),90,1);
z0star = squeeze([zp25(1,1,:),zp50(1,1,:),zp75(1,1,:),zp90(1,1,:)])';
h0 = [x0star(1,:);x0star(2,:);x0star(3,:);x0star(6,:);x0star(7,:);z0star(:,2)';t0]';

b_sape_x1 = zeros([K(2),T,inter]);
b_sape_x2 = zeros([K(2),T,inter]);
b_sape_t1 = zeros([1,T,inter]);
b_sape_t2 = zeros([1,T,inter]);
b_sape_h = zeros([K(3),T,inter]);

for i = 1:inter
beta1hat = b_thetahat(1:K(2),i);
beta2hat = b_thetahat(K(2)+1:2*K(2),i);
pihat = b_thetahat(2*K(2)+1:2*K(2)+K(3),i);
alpha1hat = b_thetahat(length(thetahat)-5,i);
alpha2hat = b_thetahat(length(thetahat)-4,i);
lambda1hat = b_thetahat(length(thetahat)-3,i);
lambda2hat = b_thetahat(length(thetahat)-2,i);
zeta1hat = b_thetahat(length(thetahat)-1,i);
zeta2hat = b_thetahat(length(thetahat),i);

dit1_0 = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + (lambda1hat*h0*pihat)';
dit2_0 = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + (lambda2hat*h0*pihat)';

sims = 1000;
u2 = normrnd(0,1,[1,1,sims]);

Lprime1_0 = (exp(dit1_0 + zeta1hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit1_0+zeta1hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;
Lprime2_0 = (exp(dit2_0 + zeta2hat*u2).*(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)) - exp(2*(dit2_0+zeta2hat*u2)))./(1 + exp(dit1_0+zeta1hat*u2) + exp(dit2_0 + zeta2hat*u2)).^2;

%x's
b_sape_x1(:,:,i) = beta1hat.*mean(Lprime1_0,3);
b_sape_x2(:,:,i) = beta2hat.*mean(Lprime2_0,3);

%h
delLdelht = lambda1hat*reshape(pihat,1,1,[]).*mean(Lprime1_0,3);
b_sape_h(:,:,i) = reshape(mean(delLdelht,1),4,13)';

%t's
ditt11 = sum(x0star.*beta1hat,1)  + alpha1hat + (lambda1hat*h0*pihat)';
ditt12 = sum(x0star.*beta2hat,1)  + alpha2hat + (lambda2hat*h0*pihat)';

ditt01 = sum(x0star.*rbeta1hat,1)  + (lambda1hat*h0*pihat)';
ditt02 = sum(x0star.*beta2hat,1)  + (lambda2hat*h0*pihat)';

L1t0 = exp(ditt01+zeta1hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L1t1 = exp(ditt11+zeta1hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

L2t0 = exp(ditt02+zeta2hat*u2)./(1 + exp(ditt01+zeta1hat*u2) + exp(ditt02 + zeta2hat*u2));
L2t1 = exp(ditt12+zeta2hat*u2)./(1 + exp(ditt11+zeta1hat*u2) + exp(ditt12 + zeta2hat*u2));

b_sape_t1(:,:,i) = mean(mean(L1t1-L1t0,3),1);
b_sape_t2(:,:,i) = mean(mean(L2t1-L2t0,3),1);
end

%Bootstraped Standard Errors
%x's
sape_x1_sd = sqrt(sum((b_sape_x1 - mean(b_sape_x1)).^2,3)./(inter-1));
sape_x2_sd = sqrt(sum((b_sape_x2 - mean(b_sape_x2)).^2,3)./(inter-1));

%t's
sape_t1_sd = sqrt(sum((b_sape_t1 - mean(b_sape_t1)).^2,3)./(inter-1));
sape_t2_sd = sqrt(sum((b_sape_t2 - mean(b_sape_t2)).^2,3)./(inter-1));

%h
sape_h_sd = sqrt(sum((b_sape_h - mean(b_sape_h)).^2,3)./(inter-1));

%Bootstrap with cre but no adjustment for endogeneity*****
%Setting things up
inter = 50;
samps = 100;
b_thetahat = zeros([size(thetahat,1),inter]);



%% Bootstrap Theta

rng(1)
s_index = randi([1,size(x,1)],samps,1,inter);
for i = 1:inter
    
    x_samp = x(s_index(:,:,i),:);
    x_samp = reshape(x_samp,samps,T,[]);
    t_samp = t(s_index(:,:,i),:);
    h_samp = h(s_index(:,:,i),:);
    fracres_temp = fracres(s_index(:,:,i),:,:);
    
    
    initial = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
    W = eye(size(initial,1));
   
    [tempthetahat] = fminsearch(@(theta) SGMM_Empiric_S2_no_endog(theta,x_samp,h_samp,t_samp,fracres_temp,W,K),initial,options);
    b_thetahat(:,i) = tempthetahat;
end

%% Bootstrap APE
b_ape_x1 = zeros([K(1),T,inter]);
b_ape_x2 = zeros([K(1),T,inter]);
b_ape_t1 = zeros([1,T,inter]);
b_ape_t2 = zeros([1,T,inter]);

for i = 1:inter
h_samp = h(s_index(:,:,i),:); 
beta1hat = b_thetahat(1:K(1),i);
beta2hat = b_thetahat(K(1)+1:2*K(1),i);
pihat = b_thetahat(2*K(1)+1:2*K(1)+K(2),i);
alpha1hat = b_thetahat(length(thetahat)-3,i);
alpha2hat = b_thetahat(length(thetahat)-2,i);
lambda1hat = b_thetahat(length(thetahat)-1,i);
lambda2hat = b_thetahat(length(thetahat),i);

idelLdelx1c_draws = ones([K(1),4,samps]);
idelLdelx2c_draws = ones([K(1),4,samps]);
iSigmal1c_mean_draws = ones([1,4,samps]);
iSigmal2c_mean_draws = ones([1,4,samps]);

rng(1) %set seed


for j=1:samps    
dit1_c = sum(x0star.*reshape(beta1hat,1,1,[]),3) + alpha1hat.*t0 + lambda1hat*h_samp(j,:)*pihat;
dit2_c = sum(x0star.*reshape(beta2hat,1,1,[]),3) + alpha2hat.*t0 + lambda2hat*h_samp(j,:)*pihat;

Lprime1_c = (exp(dit1_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit1_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
Lprime2_c = (exp(dit2_c).*(1 + exp(dit1_c) + exp(dit2_c)) - exp(2*(dit2_c)))./(1 + exp(dit1_c) + exp(dit2_c)).^2;
delLdelx1c = beta1hat.*Lprime1_c; 
idelLdelx1c_draws(:,:,j) = delLdelx1c;
delLdelx2c = beta2hat.*Lprime2_c; 
idelLdelx2c_draws(:,:,j) = delLdelx2c;

ditt11c = sum(x0star.*beta1hat,1)  + alpha1hat + lambda1hat*h_samp(j,:)*pihat;
ditt12c = sum(x0star.*beta2hat,1)  + alpha2hat + lambda2hat*h_samp(j,:)*pihat;

ditt01c = sum(x0star.*beta1hat,1)  + lambda1hat*h_samp(j,:)*pihat;
ditt02c = sum(x0star.*beta2hat,1)  + lambda2hat*h_samp(j,:)*pihat;

L1t0c = exp(ditt01c)./(1 + exp(ditt01c) + exp(ditt02c));
L1t1c = exp(ditt11c)./(1 + exp(ditt11c) + exp(ditt12c));

L2t0c = exp(ditt02c)./(1 + exp(ditt01c) + exp(ditt02c));
L2t1c = exp(ditt12c)./(1 + exp(ditt11c) + exp(ditt12c));

Sigmal1_meanc = L1t1c-L1t0c;
iSigmal1c_mean_draws(:,:,j) = Sigmal1_meanc;
Sigmal2_meanc = L2t1c-L2t0c;
iSigmal2c_mean_draws(:,:,j) = Sigmal2_meanc;
end

b_ape_x1(:,:,i) = mean(idelLdelx1c_draws,3);
b_ape_x2(:,:,i) = mean(idelLdelx2c_draws,3);
b_ape_t1(:,:,i) = mean(iSigmal1c_mean_draws,3);
b_ape_t2(:,:,i) = mean(iSigmal2c_mean_draws,3);
end

%Bootstrap Standard Errors

%x's
ape_x1_sd = sqrt(sum((b_ape_x1 - mean(b_ape_x1)).^2,3)./(inter-1));
ape_x2_sd = sqrt(sum((b_ape_x2 - mean(b_ape_x2)).^2,3)./(inter-1));

%t's
ape_t1_sd = sqrt(sum((b_ape_t1 - mean(b_ape_t1)).^2,3)./(inter-1));
ape_t2_sd = sqrt(sum((b_ape_t2 - mean(b_ape_t2)).^2,3)./(inter-1));
