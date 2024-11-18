function [k0,fk,RelError] = OptimK(k0,fk)

% OptimK
%   Calculer la valeur optimale de k0 avec la méthode des moindres carrées 
%
% Syntaxe
%   [k0,fk,RelError] = OptimK(k0,fk);
%
% Description
%   k0 : valeur initiale du vecteur d'onde, réactialisée 
%   fk : valeurs de champ pour 3 valeurs de k0, réactulisée
%
%   RelError : erreur relative sur k0

% Date de la dernière version : 25 janvier 2019
% Auteur : Mondher Besbes (LCF / CNRS / IOGS)

Mat = [k0(:).*fk(:) fk(:) ones(3,1)];
%
if rcond(Mat) < 10*eps
    ABk = InvMat(Mat)*k0(:);
else
    ABk = Mat\k0(:);
end
%
P = find(abs(fk)<max(abs(fk)));
[~,I] = sort(abs(fk(P)));
[fk(1:2),k0] = deal(fk(P(I)),[k0(P(I)) ABk(3)]);

RelError = max(abs(diff(abs([k0 k0(1)]))))/max(abs(k0));

% P = abs(fk)<max(abs(fk));
% [fk(1:2),k0] = deal(fk(P),[k0(P) ABk(3)]);
% RelError = max(abs(diff(abs(k0))))/max(abs(k0));

end