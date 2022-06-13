function [fval, psi_delta] = GMM_delta_reduced(delta,z,x,t,W)

delz = delta(1);
delx = delta(2);

q = z*delz + x*delx;

%-----delz--------
Sdelzt = zeros(size(x));
for p = 1:size(x,2)
    Sdelzt(:,p) = ((z(:,p).* normpdf(q(:,p))).*(t(:,p)-normcdf(q(:,p))))./(normcdf(q(:,p)).*(1-normcdf(q(:,p))));
end
Sdelz = sum(Sdelzt,2);

%-----delx--------
Sdelxt = zeros(size(x));
for p = 1:size(x,2)
    Sdelxt(:,p) = ((x(:,p).* normpdf(q(:,p))).*(t(:,p)-normcdf(q(:,p))))./(normcdf(q(:,p)).*(1-normcdf(q(:,p))));
end
Sdelx = sum(Sdelxt,2);



psi_delta = [Sdelz';Sdelx'];
mpsi = [mean(Sdelz);mean(Sdelx)];

fval = mpsi'*W*mpsi;
end