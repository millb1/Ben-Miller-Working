function [fval, psi_delta] = SGMM_Empiric_S1(delta,z,x,K,t,W)

%z instruments
%x varibales, dimensions: NxTxK 


delz = delta(1:K);
delx = delta(K+1:end);

q = sum(z.*reshape(delz,1,1,[]),3) + sum(x.*reshape(delx,1,1,[]),3);

%-----delz--------
Sdelzt = zeros(size(z));
for p = 1:size(z,2)
    Sdelzt(:,p) = ((z(:,p).* normpdf(q(:,p))).*(t(:,p)-normcdf(q(:,p))))./(normcdf(q(:,p)).*(1-normcdf(q(:,p))));
end
Sdelz = sum(Sdelzt,2);

%-----delx--------
Sdelxt = zeros(size(x));
for i = 1:size(x,3)
    for p = 1:size(x,2)
        Sdelxt(:,p,i) = ((x(:,p,i).* normpdf(q(:,p))).*(t(:,p)-normcdf(q(:,p))))./(normcdf(q(:,p)).*(1-normcdf(q(:,p))));
    end
end

Sdelz = sum(Sdelzt,2);
Sdelx = sum(Sdelxt,2);

psi_delta = [reshape(Sdelz,[size(z,1),size(z,3)]),reshape(Sdelx,[size(x,1),size(x,3)])]';
mpsi = mean(psi_delta,2);

fval = mpsi'*W*mpsi;
end

