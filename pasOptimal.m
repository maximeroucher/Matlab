function [alphaOptimal,fval]=pasOptimal(r0,X)
%============================================================
%
% On minimise par raaport à r la fonction :
%       \psi(r) = J(X - r \nabla J(X)).
%
%       La methode de fiminsearch est utilisée
%============================================================
x = X(1);
y = X(2);
%
% Calcul du gradient de J en X
% ----------------------------
%[dJx,dJy]=GradJ(X);
%
%psi = @(r) -2*(x-r*dJx).*(y-r*dJy)-4*(x-r*dJx)+(x-r*dJx).^2+2*(y-r*dJy).^2 ;
%
%[alphaOptimal,fval] = fminsearch(psi,r0);
%
% Calcul du gradient de J en X
% ----------------------------
[dJx,dJy]=GradJ(X);
%
% Cas de la fonction J
% --------------------
phi = @(r)((X(1)-r.*dJx)-1).^2+10.*((X(1)-r.*dJx).^2-(X(2)-r.*dJy)).^2;
%
[alphaOptimal,fval] = fminsearch(phi,r0);
%
end

