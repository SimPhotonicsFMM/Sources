function AddIndex(varargin)

% 
% AddIndex 
%   Définir un nouveau matériau
%
% Syntaxe
%   AddIndex('Param1',Val1,'Param2',Val2,...);
%
% Description
%   Param     : Nom du champ
%   Val       : Valeur associée à Param

% Date de la dernière version : 03 mai 2021
% Auteur : Mondher Besbes (LCF / CNRS / IOGS)

load('IndexData.mat','Index')   % Index structure

f = varargin(1:2:end); % nom des champs à modifier
c = varargin(2:2:end); % valeurs associées

FieldData = fieldnames(Index(1));


TabName = {Index.Name};
Pn = find(ismember(f,'Name'));
P = find(ismember(TabName,c(Pn)));
if isempty(P)
    N = length(Index)+1; % Nouveau matériau
else
    N = P;
end
%
for k = 1:length(f)
    P = ismember(FieldData,f(k));
    if any(P) 
        Index(N).(FieldData{P}) = c{k};
    else % nouveau champ
        Index(N).(f{k}) = c{k};
        warning(['Nouvelle donnée : ', f{k}]); 
    end
end

% Ranger par ordre alphabétique
TabName = {Index.Name}; 
[~,P] = sort(TabName);
Index = Index(P);

save IndexData.mat Index

end

