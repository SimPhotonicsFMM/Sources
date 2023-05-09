function [a,b] = CoefPade(x,y,n,m)

% CoefPade 
%   Calcul des coefficients a et b des polynômes de Padé pour des couples
%   de points (xi,yi)
%
%           Pn,m(x) = (a0+a1*x+...+an*x^n)/(1+b1*x+...+bn*x^m);
%
% Syntaxe
%   [a,b] = CoefPade(x,y,n,m);
%
% Description
%   x,y : valeurs de la fonction échantiollonée
%   n,m : ordres des polynômes de Padé
%
%   a,b : coefficients des polynômes de Padé

% Date de la dernière version : 26 septembre 2013 
% Auteur : Mondher Besbes (LCFIO / CNRS)

A = ones(length(y),n+m+1);

for k = 1:n,
    A(:,k+1) = x(:).^k;
end
%
for k = 1:m,
    A(:,k+n+1) = -y(:).*x(:).^k;
end

ab = A\y(:);
%ab = pinv(A)*y(:);
a = ab(1:n+1);
b = ab(n+2:end);


return