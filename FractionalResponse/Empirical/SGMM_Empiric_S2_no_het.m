function[fval, psi_theta] = SGMM_Empiric_S2_no_het(theta,x,t,fracres,W,K)

 %K - number of variables in x


    
%Naming coefficients
beta1 = theta(1:K);
beta2 = theta(K+1:2*K);
alpha1 = theta(length(theta)-1);
alpha2 = theta(length(theta));





dit1 = sum(x.*reshape(beta1,1,1,[]),3) + t*alpha1 ;
dit2 = sum(x.*reshape(beta2,1,1,[]),3) + t*alpha2 ;


L1 = exp(dit1)./(1 + exp(dit1) + exp(dit2));
L2 = exp(dit2)./(1 + exp(dit1) + exp(dit2));

Lprime1 = (exp(dit1).*(1 + exp(dit1) + exp(dit2)) - exp(2*(dit1)))./(1 + exp(dit1) + exp(dit2)).^2;
Lprime2 = (exp(dit2).*(1 + exp(dit1) + exp(dit2)) - exp(2*(dit2)))./(1 + exp(dit1) + exp(dit2)).^2;

M1 = (-exp(dit1).*exp(dit2))./((1 + exp(dit1) + exp(dit2)).^2);
M2 = (-exp(dit2).*exp(dit1))./((1 + exp(dit1) + exp(dit2)).^2);

P1 = (-exp(dit1))./((1 + exp(dit1) + exp(dit2)).^2);
P2 = (-exp(dit2))./((1 + exp(dit1) + exp(dit2)).^2);

oneminussum = 1- (L1+L2);



%Calculation of integrals
int_oneminussum_inf = oneminussum ;
oneminussumintegral_infty = int_oneminussum_inf;



int_L1_inf = L1 ;
L1integral_infty = int_L1_inf;

int_L2_inf = L2 ;
L2integral_infty = int_L2_inf;



%Lprime integrals with x's
int_Lprime1_inf_x = Lprime1.*reshape(x,[size(x,1),size(x,2),1,size(x,3)]);
Lprime1integral_infty_x = int_Lprime1_inf_x ;

int_Lprime2_inf_x = Lprime2.*reshape(x,[size(x,1),size(x,2),1,size(x,3)]);
Lprime2integral_infty_x = int_Lprime2_inf_x ;



%M integrals with x's
int_M1_inf_x = M1 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) ;
M1integral_infty_x = int_M1_inf_x ;

int_M2_inf_x = M2 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) ;
M2integral_infty_x = int_M2_inf_x;


%P integral with x's
int_P1_inf_x = P1 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) ;
P1integral_infty_x = int_P1_inf_x ;

int_P2_inf_x = P2 .* reshape(x,[size(x,1),size(x,2),1,size(x,3)]) ;
P2integral_infty_x = int_P2_inf_x ;



%t integrals
int_Lprime1_inf_t = Lprime1 .* t ;
Lprime1integral_infty_t = int_Lprime1_inf_t ;

int_Lprime2_inf_t = Lprime2 .* t ;
Lprime2integral_infty_t = int_Lprime2_inf_t ;



int_M1_inf_t = M1 .* t ;
M1integral_infty_t = int_M1_inf_t ;

int_M2_inf_t = M2 .* t ;
M2integral_infty_t = int_M2_inf_t ;



int_P1_inf_t = P1 .* t ;
P1integral_infty_t = int_P1_inf_t ;

int_P2_inf_t = P2 .* t ;
P2integral_infty_t = int_P2_inf_t ;






%beta1

Sbeta1t = zeros(size(x));
for i = 1:size(x,3)
    for p = 1:size(x,2)
    Sbeta1t(:,p,i) =fracres(:,p,2).*(M2integral_infty_x(:,p,i)./L2integral_infty(:,p)) ...
               + fracres(:,p,1).*(Lprime1integral_infty_x(:,p,i)./L1integral_infty(:,p)) ...
               + fracres(:,p,3).*(P1integral_infty_x(:,p,i)./(oneminussumintegral_infty(:,p)));
    end
end    
Sbeta1 = reshape(sum(Sbeta1t,2),size(x,1),[]);

%beta 2
Sbeta2t = zeros(size(x));
for i = 1:size(x,3)
    for p = 1:size(x,2)
    Sbeta2t(:,p,i) = fracres(:,p,1).*(M1integral_infty_x(:,p,i)./L1integral_infty(:,p)) ...
               + fracres(:,p,2).*(Lprime2integral_infty_x(:,p,i)./L2integral_infty(:,p)) ...
               + fracres(:,p,3).*(P2integral_infty_x(:,p,i)./(oneminussumintegral_infty(:,p)));
    end
end
Sbeta2 = reshape(sum(Sbeta2t,2),size(x,1),[]);

%-----alpha 1--------
Salpha1t = zeros(size(x,1),size(x,2));
for p = 1:size(x,2)
    Salpha1t(:,p) = fracres(:,p,2).*(M2integral_infty_t(:,p)./L2integral_infty(:,p)) ...
               + fracres(:,p,1).*(Lprime1integral_infty_t(:,p)./L1integral_infty(:,p)) ...
               + fracres(:,p,3).*(P1integral_infty_t(:,p)./(oneminussumintegral_infty(:,p)));
end
Salpha1 = sum(Salpha1t,2);

%-----alpha 2--------
Salpha2t = zeros(size(x,1),size(x,2));
for p = 1:size(x,2)
    Salpha2t(:,p) = fracres(:,p,1).*(M1integral_infty_t(:,p)./L1integral_infty(:,p)) ...
               + fracres(:,p,2).*(Lprime2integral_infty_t(:,p)./L2integral_infty(:,p)) ...
               + fracres(:,p,3).*(P2integral_infty_t(:,p)./(oneminussumintegral_infty(:,p)));
end
Salpha2 = sum(Salpha2t,2);



psi_theta = [reshape(Sbeta1,size(x,1),[]),reshape(Sbeta2,size(x,1),[]),Salpha1,Salpha2]';
mpsi = mean(psi_theta');
fval = mpsi*W*mpsi';
end
