function[fval, psi_theta] = SGMM_Empiric_S2(theta,delta,z,x,h,t,u,w,fracres,W,pie,K)

%K-vector of variable lengths
    %K(1) - number of varibales in z
    %K(2) - number of variables in x
    %K(3) - number of variables in h

    
%Naming coefficients
delz = delta(1:K(1));
delx = delta(K(1)+1:end);
beta1 = theta(1:K(2));
beta2 = theta(K(2)+1:2*K(2));
pi = theta(2*K(2)+1:2*K(2)+K(3));
alpha1 = theta(length(theta)-5);
alpha2 = theta(length(theta)-4);
lambda1 = theta(length(theta)-3);
lambda2 = theta(length(theta)-2);
zeta1 = theta(length(theta)-1);
zeta2 = theta(length(theta));

%
q = sum(z.*reshape(delz,1,1,[]),3) + sum(x.*reshape(delx,1,1,[]),3);

dit1 = sum(x.*reshape(beta1,1,1,[]),3) + t*alpha1 + lambda1*h*pi;
dit2 = sum(x.*reshape(beta2,1,1,[]),3) + t*alpha2 + lambda2*h*pi;

% For integrals from -q to infinity
L1 = (exp(u).*exp(zeta1.*((dit1./zeta1)+u-q-(((u-q).^2)./(2.*zeta1))))./(sqrt(2.*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q)))));
L2 = (exp(u).*exp(zeta2.*((dit2./zeta2)+u-q-(((u-q).^2)./(2.*zeta2)))))./(sqrt(2.*pie)*(1 + exp(zeta1*((dit1/zeta1)+u-q)) + exp(zeta2*((dit2/zeta2)+u-q))));

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

%Lprime integrals with x's
int_Lprime1_inf_x = Lprime1.*reshape(x,[size(x,1),size(x,2),1,size(x,3)]).* w;
Lprime1integral_infty_x = reshape(sum(int_Lprime1_inf_x,3),size(x));

int_Lprime2_inf_x = Lprime2.*reshape(x,[size(x,1),size(x,2),1,size(x,3)]).* w;
Lprime2integral_infty_x = reshape(sum(int_Lprime2_inf_x,3),size(x));

int_Lprime1_neginf_x = Lprime1_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
Lprime1integral_neginfty_x = reshape(sum(int_Lprime1_neginf_x,3),size(x));

int_Lprime2_neginf_x = Lprime2_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
Lprime2integral_neginfty_x = reshape(sum(int_Lprime2_neginf_x,3),size(x));

%M integrals with x's
int_M1_inf_x = M1 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
M1integral_infty_x = reshape(sum(int_M1_inf_x,3),size(x));

int_M2_inf_x = M2 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
M2integral_infty_x = reshape(sum(int_M2_inf_x,3),size(x));

int_M1_neginf_x = M1_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
M1integral_neginfty_x = reshape(sum(int_M1_neginf_x,3),size(x));

int_M2_neginf_x = M2_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
M2integral_neginfty_x = reshape(sum(int_M2_neginf_x,3),size(x));

%P integral with x's
int_P1_inf_x = P1 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
P1integral_infty_x = reshape(sum(int_P1_inf_x,3),size(x));

int_P2_inf_x = P2 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
P2integral_infty_x = reshape(sum(int_P2_inf_x,3),size(x));

int_P1_neginf_x = P1_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
P1integral_neginfty_x = reshape(sum(int_P1_neginf_x,3),size(x));

int_P2_neginf_x = P2_neginf .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) .* w;
P2integral_neginfty_x = reshape(sum(int_P2_neginf_x,3),size(x));

%t integrals
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

%hpi integrals

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

%u-q integrals

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

%lambda integrals

int_expo_inf_lambdah = expo .* reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
expointegral_infty_lambdah = reshape(sum(int_expo_inf_lambdah,3),[size(expo,1),size(expo,2),size(h,2)]);

int_expo_neginf_lambdah = expo_neginf .* reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
expointegral_neginfty_lambdah = sum(int_expo_neginf_lambdah,3);

int_P1_inf_lambdah = -P1 .* lambda1.*reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
P1integral_infty_lambdah = reshape(sum(int_P1_inf_lambdah,3),[size(expo,1),size(expo,2),size(h,2)]);

int_P1_neginf_lambdah = -P1_neginf .* lambda1.*reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
P1integral_neginfty_lambdah = reshape(sum(int_P1_neginf_lambdah,3),[size(expo,1),size(expo,2),size(h,2)]);

int_P2_inf_lambdah = -P2 .* lambda2.*reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
P2integral_infty_lambdah = reshape(sum(int_P2_inf_lambdah,3),[size(expo,1),size(expo,2),size(h,2)]);

int_P2_neginf_lambdah = -P2_neginf .* lambda2.*reshape(h,[size(h,1),1,size(h,3),size(h,2)]) .* w;
P2integral_neginfty_lambdah = reshape(sum(int_P2_neginf_lambdah,3),[size(expo,1),size(expo,2),size(h,2)]);

%beta1

Sbeta1t = zeros(size(x));
for i = 1:size(x,3)
    for p = 1:size(x,2)
    Sbeta1t(:,p,i) =t(:,p).*fracres(:,p,2).*(M2integral_infty_x(:,p,i)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,1).*(Lprime1integral_infty_x(:,p,i)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P1integral_infty_x(:,p,i)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,2).*(M2integral_neginfty_x(:,p,i)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,1).*(Lprime1integral_neginfty_x(:,p,i)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P1integral_neginfty_x(:,p,i)./(oneminussumintegral_neginfty(:,p)));
    end
end    
Sbeta1 = reshape(sum(Sbeta1t,2),size(x,1),[]);

%beta 2
Sbeta2t = zeros(size(x));
for i = 1:size(x,3)
    for p = 1:size(x,2)
    Sbeta2t(:,p,i) = t(:,p).*fracres(:,p,1).*(M1integral_infty_x(:,p,i)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_x(:,p,i)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_x(:,p,i)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_x(:,p,i)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_x(:,p,i)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_x(:,p,i)./(oneminussumintegral_neginfty(:,p)));
    end
end
Sbeta2 = reshape(sum(Sbeta2t,2),size(x,1),[]);

%-----alpha 1--------
Salpha1t = zeros(size(x,1),size(x,2));
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
Salpha2t = zeros(size(x,1),size(x,2));
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
Slambda1t = zeros(size(x,1),size(x,2));
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
Slambda2t = zeros(size(x,1),size(x,2));
for p = 1:size(x,2)
    Slambda2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_pih(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_pih(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_pih(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_pih(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_pih(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_pih(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Slambda2 = sum(Slambda2t,2);

%pi
Spi1t = zeros(size(x,1),size(x,2),size(h,2));
for i = 1:K(3)
    for p = 1:size(x,2)
Spi1t(:,p,i) = t(:,p).*fracres(:,p,1).*(P1integral_infty_lambdah(:,p,i)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(P2integral_infty_lambdah(:,p,i)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(expointegral_infty_lambdah(:,p,i)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(P1integral_neginfty_lambdah(:,p,i)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(P2integral_neginfty_lambdah(:,p,i)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(expointegral_neginfty_lambdah(:,p,i)./(oneminussumintegral_neginfty(:,p)));
    end
end
Spi = reshape(sum(Spi1t,2),size(x,1),[]);

%-----zeta 1--------
Szeta1t = zeros(size(x,1),size(x,2));
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
Szeta2t = zeros(size(x,1),size(x,2));
for p = 1:size(x,2)
    Szeta2t(:,p) = t(:,p).*fracres(:,p,1).*(M1integral_infty_u(:,p)./L1integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,2).*(Lprime2integral_infty_u(:,p)./L2integral_infty(:,p)) ...
               + t(:,p).*fracres(:,p,3).*(P2integral_infty_u(:,p)./(oneminussumintegral_infty(:,p))) ...
               + (1-t(:,p)).*fracres(:,p,1).*(M1integral_neginfty_u(:,p)./L1integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,2).*(Lprime2integral_neginfty_u(:,p)./L2integral_neginfty(:,p)) ...
               + (1-t(:,p)).*fracres(:,p,3).*(P2integral_neginfty_u(:,p)./(oneminussumintegral_neginfty(:,p)));
end
Szeta2 = sum(Szeta2t,2);


psi_theta = [reshape(Sbeta1,size(x,1),[]),reshape(Sbeta2,size(x,1),[]),reshape(Spi,size(x,1),[]),Salpha1,Salpha2,Slambda1,Slambda2,Szeta1,Szeta2]';
mpsi = mean(psi_theta');
fval = mpsi*W*mpsi';
end

