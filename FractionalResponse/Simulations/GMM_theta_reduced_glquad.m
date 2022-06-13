function [fval,psi_theta] = GMM_theta_reduced_glquad(theta,delta,z,x,h,u,t,w,fracres,W,pie)
delz = delta(1);
delx = delta(2);
beta1 = theta(1);
beta2 = theta(2);
alpha1 = theta(3);
alpha2 = theta(4);
lambda1 = theta(5);
lambda2 = theta(6);
pi = theta(7:9)';
zeta1 = theta(10);
zeta2 = theta(11);

q = z*delz + x*delx;

dit1 = x*beta1+t*alpha1+lambda1*h*pi;
dit2 = x*beta2+t*alpha2+lambda2*h*pi;
% For integrals from -q to infinity
L1 = (exp(u).*exp(zeta1.*((dit1./zeta1)+u-q-(((u-q).^2)./(2.*zeta1)))))./(sqrt(2.*pie).*(1 + exp(zeta1.*((dit1./zeta1)+u-q)) + exp(zeta2.*((dit2./zeta2)+u-q))));
L2 = (exp(u).*exp(zeta2*((dit2/zeta2)+u-q-(((u-q).^2)/(2*zeta2)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))));

Lprime1 = (exp(u).*(exp((zeta1+zeta2)*(((dit1+dit2)/(zeta1+zeta2))+u-q-(((u-q).^2)/(2*(zeta1+zeta2)))))+exp(zeta1*((dit1/zeta1)+u-q-(((u-q).^2)/(2*zeta1))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));
Lprime2 = (exp(u).*(exp((zeta1+zeta2)*(((dit1+dit2)/(zeta1+zeta2))+u-q-(((u-q).^2)/(2*(zeta1+zeta2)))))+exp(zeta2*((dit2/zeta2)+u-q-(((u-q).^2)/(2*zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));

M1 = -(exp(u).*exp((zeta1+zeta2)*(((dit1+dit2)./(zeta1+zeta2))+u-q-(((u-q).^2)/(2*(zeta1+zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));
M2 = -(exp(u).*exp((zeta1+zeta2)*(((dit1+dit2)./(zeta1+zeta2))+u-q-(((u-q).^2)/(2*(zeta1+zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));

P1 = -(exp(u).*exp(zeta1*((dit1/zeta1)+u-q-(((u-q).^2)/(2*zeta1)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));
P2 = -(exp(u).*exp(zeta2*((dit2/zeta2)+u-q-(((u-q).^2)/(2*zeta2)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));

expo = -(exp(u).*(lambda1*exp(zeta1*((dit1/zeta1)+u-q-(((u-q).^2)/(2*zeta1))))+lambda2*exp(zeta2*((dit2/zeta2)+u-q-(((u-q).^2)/(2*zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))).^(2));
oneminussum = (exp(u).*exp(-(((u-q).^2)/(2))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))));
%For integrals from -infinity to -q
L1_neginf = (exp(u).*exp(zeta1*((dit1/zeta1)-u-q-(((-u-q).^2)/(2*zeta1)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))));
L2_neginf = (exp(u).*exp(zeta2*((dit2/zeta2)-u-q-(((-u-q).^2)/(2*zeta2)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))));

Lprime1_neginf = (exp(u).*(exp((zeta1+zeta2)*(((dit1+dit2)/(zeta1+zeta2))-u-q-(((-u-q).^2)/(2*(zeta1+zeta2)))))+exp(zeta1*((dit1/zeta1)-u-q-(((-u-q).^2)/(2*zeta1))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));
Lprime2_neginf = (exp(u).*(exp((zeta1+zeta2)*(((dit1+dit2)/(zeta1+zeta2))-u-q-(((-u-q).^2)/(2*(zeta1+zeta2)))))+exp(zeta2*((dit2/zeta2)-u-q-(((-u-q).^2)/(2*zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));

M1_neginf = -(exp(u).*exp((zeta1+zeta2)*(((dit1+dit2)./(zeta1+zeta2))-u-q-(((-u-q).^2)/(2*(zeta1+zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));
M2_neginf = -(exp(u).*exp((zeta1+zeta2)*(((dit1+dit2)./(zeta1+zeta2))-u-q-(((-u-q).^2)/(2*(zeta1+zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));

P1_neginf = -(exp(u).*exp(zeta1*((dit1/zeta1)-u-q-(((-u-q).^2)/(2*zeta1)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));
P2_neginf = -(exp(u).*exp(zeta2*((dit2/zeta2)-u-q-(((-u-q).^2)/(2*zeta2)))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));

expo_neginf = -(exp(u).*(lambda1*exp(zeta1*((dit1/zeta1)-u-q-(((-u-q).^2)/(2*zeta1))))+lambda2*exp(zeta2*((dit2/zeta2)-u-q-(((-u-q).^2)/(2*zeta2))))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))).^(2));

oneminussum_neginf = (exp(u).*exp(-(((-u-q).^2)/(2))))./(sqrt(2*pie)*(1 + exp(zeta1*((dit1/zeta1)-u-q)) + exp(zeta2*((dit2/zeta2)-u-q))));

%qplus = double((u >= -q)); 
%qminus = double((u < -q));  

%Calculation of integrals


int_oneminussum_inf = oneminussum .* w;
oneminussumintegral_infty = sum(int_oneminussum_inf,3);

int_oneminussum_neginf = oneminussum_neginf .* w;
oneminussumintegral_neginfty = sum(int_oneminussum_neginf,3);

int_L1_inf = L1 .* w;
L1integral_infty = sum(int_L1_inf,3);

int_L2_inf = L2 .* w;
L2integral_infty = sum(int_L2_inf,3);

int_L1_neginf = L1_neginf .* w;
L1integral_neginfty = sum(int_L1_neginf,3);


int_L2_neginf = L2_neginf .* w;
L2integral_neginfty = sum(int_L2_neginf,3);

int_Lprime1_inf_x = Lprime1 .* x .* w;
Lprime1integral_infty_x = sum(int_Lprime1_inf_x,3);

int_Lprime2_inf_x = Lprime2 .* x .* w;
Lprime2integral_infty_x = sum(int_Lprime2_inf_x,3);

int_Lprime1_neginf_x = Lprime1_neginf .* x .* w;
Lprime1integral_neginfty_x = sum(int_Lprime1_neginf_x,3);

int_Lprime2_neginf_x = Lprime2_neginf .* x .* w;
Lprime2integral_neginfty_x = sum(int_Lprime2_neginf_x,3);

int_M1_inf_x = M1 .* x .* w;
M1integral_infty_x = sum(int_M1_inf_x,3);

int_M2_inf_x = M2 .* x .* w;
M2integral_infty_x = sum(int_M2_inf_x,3);

int_M1_neginf_x = M1_neginf .* x .* w;
M1integral_neginfty_x = sum(int_M1_neginf_x,3);

int_M2_neginf_x = M2_neginf .* x .* w;
M2integral_neginfty_x = sum(int_M2_neginf_x,3);

int_P1_inf_x = P1 .* x .* w;
P1integral_infty_x = sum(int_P1_inf_x,3);

int_P2_inf_x = P2 .* x .* w;
P2integral_infty_x = sum(int_P2_inf_x,3);

int_P1_neginf_x = P1_neginf .* x .* w;
P1integral_neginfty_x = sum(int_P1_neginf_x,3);

int_P2_neginf_x = P2_neginf .* x .* w;
P2integral_neginfty_x = sum(int_P2_neginf_x,3);

int_Lprime1_inf_t = Lprime1 .* t .* w;
Lprime1integral_infty_t = sum(int_Lprime1_inf_t,3);

int_Lprime2_inf_t = Lprime2 .* t .* w;
Lprime2integral_infty_t = sum(int_Lprime2_inf_t,3);

int_Lprime1_neginf_t = Lprime1_neginf .* t .* w;
Lprime1integral_neginfty_t = sum(int_Lprime1_neginf_t,3);

int_Lprime2_neginf_t = Lprime2_neginf .* t .* w;
Lprime2integral_neginfty_t = sum(int_Lprime2_neginf_t,3);

int_M1_inf_t = M1 .* t .* w;
M1integral_infty_t = sum(int_M1_inf_t,3);

int_M2_inf_t = M2 .* t .* w;
M2integral_infty_t = sum(int_M2_inf_t,3);

int_M1_neginf_t = M1_neginf .* t .* w;
M1integral_neginfty_t = sum(int_M1_neginf_t,3);

int_M2_neginf_t = M2_neginf .* t .* w;
M2integral_neginfty_t = sum(int_M2_neginf_t,3);

int_P1_inf_t = P1 .* t .* w;
P1integral_infty_t = sum(int_P1_inf_t,3);

int_P2_inf_t = P2 .* t .* w;
P2integral_infty_t = sum(int_P2_inf_t,3);

int_P1_neginf_t = P1_neginf .* t .* w;
P1integral_neginfty_t = sum(int_P1_neginf_t,3);

int_P2_neginf_t = P2_neginf .* t .* w;
P2integral_neginfty_t = sum(int_P2_neginf_t,3);

int_Lprime1_inf_pih = Lprime1 .* (h*pi) .* w;
Lprime1integral_infty_pih = sum(int_Lprime1_inf_pih,3);

int_Lprime2_inf_pih = Lprime2 .* (h*pi) .* w;
Lprime2integral_infty_pih = sum(int_Lprime2_inf_pih,3);

int_Lprime1_neginf_pih = Lprime1_neginf .* (h*pi) .* w;
Lprime1integral_neginfty_pih = sum(int_Lprime1_neginf_pih,3);

int_Lprime2_neginf_pih = Lprime2_neginf .* (h*pi) .* w;
Lprime2integral_neginfty_pih = sum(int_Lprime2_neginf_pih,3);

int_M1_inf_pih = M1 .* (h*pi) .* w;
M1integral_infty_pih = sum(int_M1_inf_pih,3);

int_M2_inf_pih = M2 .* (h*pi) .* w;
M2integral_infty_pih = sum(int_M2_inf_pih,3);

int_M1_neginf_pih = M1_neginf .* (h*pi) .* w;
M1integral_neginfty_pih = sum(int_M1_neginf_pih,3);

int_M2_neginf_pih = M2_neginf .* (h*pi) .* w;
M2integral_neginfty_pih = sum(int_M2_neginf_pih,3);

int_P1_inf_pih = P1 .* (h*pi) .* w;
P1integral_infty_pih = sum(int_P1_inf_pih,3);

int_P2_inf_pih = P2 .* (h*pi) .* w;
P2integral_infty_pih = sum(int_P2_inf_pih,3);

int_P1_neginf_pih = P1_neginf .* (h*pi) .* w;
P1integral_neginfty_pih = sum(int_P1_neginf_pih,3);

int_P2_neginf_pih = P2_neginf .* (h*pi) .* w;
P2integral_neginfty_pih = sum(int_P2_neginf_pih,3);

int_Lprime1_inf_u = Lprime1 .* (u-q) .* w;
Lprime1integral_infty_u = sum(int_Lprime1_inf_u,3);

int_Lprime2_inf_u = Lprime2 .* (u-q) .* w;
Lprime2integral_infty_u = sum(int_Lprime2_inf_u,3);

int_Lprime1_neginf_u = Lprime1_neginf .* (-u-q) .* w;
Lprime1integral_neginfty_u = sum(int_Lprime1_neginf_u,3);

int_Lprime2_neginf_u = Lprime2_neginf .* (-u-q) .* w;
Lprime2integral_neginfty_u = sum(int_Lprime2_neginf_u,3);

int_M1_inf_u = M1 .* (u-q) .* w;
M1integral_infty_u = sum(int_M1_inf_u,3);

int_M2_inf_u = M2 .* (u-q) .* w;
M2integral_infty_u = sum(int_M2_inf_u,3);

int_M1_neginf_u = M1_neginf .* (-u-q) .* w;
M1integral_neginfty_u = sum(int_M1_neginf_u,3);

int_M2_neginf_u = M2_neginf .* (-u-q) .* w;
M2integral_neginfty_u = sum(int_M2_neginf_u,3);

int_P1_inf_u = P1 .* (u-q) .* w;
P1integral_infty_u = sum(int_P1_inf_u,3);

int_P2_inf_u = P2 .* (u-q) .* w;
P2integral_infty_u = sum(int_P2_inf_u,3);

int_P1_neginf_u = P1_neginf .* (-u-q) .* w;
P1integral_neginfty_u = sum(int_P1_neginf_u,3);

int_P2_neginf_u = P2_neginf .* (-u-q) .* w;
P2integral_neginfty_u = sum(int_P2_neginf_u,3);

int_expo_inf_lambdah1 = expo .* h(:,1) .* w;
expointegral_infty_lambdah1 = sum(int_expo_inf_lambdah1,3);

int_expo_inf_lambdah2 = expo .* h(:,2) .* w;
expointegral_infty_lambdah2 = sum(int_expo_inf_lambdah2,3);

int_expo_inf_lambdah4 = expo .* h(:,3) .* w;
expointegral_infty_lambdah4 = sum(int_expo_inf_lambdah4,3);

int_expo_neginf_lambdah1 = expo_neginf .* h(:,1) .* w;
expointegral_neginfty_lambdah1 = sum(int_expo_neginf_lambdah1,3);

int_expo_neginf_lambdah2 = expo_neginf .* h(:,2) .* w;
expointegral_neginfty_lambdah2 = sum(int_expo_neginf_lambdah2,3);

int_expo_neginf_lambdah4 = expo_neginf .* h(:,3) .* w;
expointegral_neginfty_lambdah4 = sum(int_expo_neginf_lambdah4,3);

int_P1_inf_lambdah1 = -P1 .* lambda1.*h(:,1) .* w;
P1integral_infty_lambdah1 = sum(int_P1_inf_lambdah1,3);

int_P1_inf_lambdah2 = -P1 .* lambda1.*h(:,2) .* w;
P1integral_infty_lambdah2 = sum(int_P1_inf_lambdah2,3);

int_P1_inf_lambdah4 = -P1 .* lambda1.*h(:,3) .* w;
P1integral_infty_lambdah4 = sum(int_P1_inf_lambdah4,3);

int_P1_neginf_lambdah1 = -P1_neginf .* lambda1.*h(:,1) .* w;
P1integral_neginfty_lambdah1 = sum(int_P1_neginf_lambdah1,3);

int_P1_neginf_lambdah2 = -P1_neginf .* lambda1.*h(:,2) .* w;
P1integral_neginfty_lambdah2 = sum(int_P1_neginf_lambdah2,3);

int_P1_neginf_lambdah4 = -P1_neginf .* lambda1.*h(:,3) .* w;
P1integral_neginfty_lambdah4 = sum(int_P1_neginf_lambdah4,3);

int_P2_inf_lambdah1 = -P2 .* lambda2.*h(:,1) .* w;
P2integral_infty_lambdah1 = sum(int_P2_inf_lambdah1,3);

int_P2_inf_lambdah2 = -P2 .* lambda2.*h(:,2) .* w;
P2integral_infty_lambdah2 = sum(int_P2_inf_lambdah2,3);

int_P2_inf_lambdah4 = -P2 .* lambda2.*h(:,3) .* w;
P2integral_infty_lambdah4 = sum(int_P2_inf_lambdah4,3);

int_P2_neginf_lambdah1 = -P2_neginf .* lambda2.*h(:,1) .* w;
P2integral_neginfty_lambdah1 = sum(int_P2_neginf_lambdah1,3);

int_P2_neginf_lambdah2 = -P2_neginf .* lambda2.*h(:,2) .* w;
P2integral_neginfty_lambdah2 = sum(int_P2_neginf_lambdah2,3);

int_P2_neginf_lambdah4 = -P2_neginf .* lambda2.*h(:,3) .* w;
P2integral_neginfty_lambdah4 = sum(int_P2_neginf_lambdah4,3);

%-----Beta 1--------
Sbeta1t = zeros(size(x));
for p = 1:size(x,2)
    Sbeta1t(:,p) =t(:,p).*fracres(:,p,2).*(M2integral_infty_x(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,1).*(Lprime1integral_infty_x(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P1integral_infty_x(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,2).*(M2integral_neginfty_x(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,1).*(Lprime1integral_neginfty_x(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P1integral_neginfty_x(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Sbeta1 = sum(Sbeta1t,2);

%-----Beta 2--------
Sbeta2t = zeros(size(x));
for p = 1:size(x,2)
    Sbeta2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_x(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_x(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_x(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_x(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_x(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_x(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Sbeta2 = sum(Sbeta2t,2);


%-----alpha 1--------
Salpha1t = zeros(size(x));
for p = 1:size(x,2)
    Salpha1t(:,p) = t(:,p).*fracres(:,p,2).*(M2integral_infty_t(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,1).*(Lprime1integral_infty_t(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P1integral_infty_t(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,2).*(M2integral_neginfty_t(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,1).*(Lprime1integral_neginfty_t(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P1integral_neginfty_t(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Salpha1 = sum(Salpha1t,2);

%-----alpha 2--------
Salpha2t = zeros(size(x));
for p = 1:size(x,2)
    Salpha2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_t(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_t(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_t(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_t(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_t(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_t(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Salpha2 = sum(Salpha2t,2);

%-----lambda 1--------
Slambda1t = zeros(size(x));
for p = 1:size(x,2)
    Slambda1t(:,p) = t(:,p).*fracres(:,p,2).*(M2integral_infty_pih(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,1).*(Lprime1integral_infty_pih(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P1integral_infty_pih(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,2).*(M2integral_neginfty_pih(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,1).*(Lprime1integral_neginfty_pih(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P1integral_neginfty_pih(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Slambda1 = sum(Slambda1t,2);

%-----lambda 2--------
Slambda2t = zeros(size(x));
for p = 1:size(x,2)
    Slambda2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_pih(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_pih(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_pih(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_pih(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_pih(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_pih(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Slambda2 = sum(Slambda2t,2);

%-----pi 1--------
Spi1t = zeros(size(x));
for p = 1:size(x,2)
Spi1t(:,p) = t(:,p).*fracres(:,p,1).*(P1integral_infty_lambdah1(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(P2integral_infty_lambdah1(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(expointegral_infty_lambdah1(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(P1integral_neginfty_lambdah1(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(P2integral_neginfty_lambdah1(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(expointegral_neginfty_lambdah1(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Spi1 = sum(Spi1t,2);

%-----pi 2--------
Spi2t = zeros(size(x));
for p = 1:size(x,2)
Spi2t(:,p) = t(:,p).*fracres(:,p,1).*(P1integral_infty_lambdah2(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(P2integral_infty_lambdah2(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(expointegral_infty_lambdah2(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(P1integral_neginfty_lambdah2(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(P2integral_neginfty_lambdah2(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(expointegral_neginfty_lambdah2(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Spi2 = sum(Spi2t,2);


%-----pi 4--------
Spi4t = zeros(size(x));
for p = 1:size(x,2)
Spi4t(:,p) = t(:,p).*fracres(:,p,1).*(P1integral_infty_lambdah4(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(P2integral_infty_lambdah4(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(expointegral_infty_lambdah4(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(P1integral_neginfty_lambdah4(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(P2integral_neginfty_lambdah4(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(expointegral_neginfty_lambdah4(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Spi4 = sum(Spi4t,2);

%-----zeta 1--------
Szeta1t = zeros(size(x));
for p = 1:size(x,2)
    Szeta1t(:,p) = t(:,p).*fracres(:,p,2).*(M2integral_infty_u(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,1).*(Lprime1integral_infty_u(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P1integral_infty_u(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,2).*(M2integral_neginfty_u(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,1).*(Lprime1integral_neginfty_u(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P1integral_neginfty_u(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Szeta1 = sum(Szeta1t,2);

%-----zeta 2--------
Szeta2t = zeros(size(x));
for p = 1:size(x,2)
    Szeta2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_u(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_u(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_u(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_u(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_u(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_u(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Szeta2 = sum(Szeta2t,2);

psi_theta = [Sbeta1';Sbeta2';Salpha1';Salpha2';Slambda1';Slambda2';Spi1';Spi2';Spi4';Szeta1';Szeta2'];

mpsi = [mean(Sbeta1);mean(Sbeta2);mean(Salpha1);mean(Salpha2);mean(Slambda1);mean(Slambda2);mean(Spi1);mean(Spi2);mean(Spi4);mean(Szeta1);mean(Szeta2)];

fval = mpsi'*W*mpsi;
end
