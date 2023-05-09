function n = IndexVal(Name,x)
 
% IndexVal 
%   Calculate the refractive index of a material from "IndexData.mat"
%
% Syntax
%   n = IndexVal(Name,x)    % List of values of n
%   n = IndexVal(Name)      % hundle function of n
%
% Description
%   Name   : Name of material
%   x     : Wavelength (µm)

% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

global Index

% Cas plusieurs matériaux
if iscell(Name)
    if nargin == 1  % handle function
        n = cell(1,length(Name));
        for k = 1:length(Name)
            n{k} = IndexVal(Name{k});
        end
    else            % Tableau de valeurs
        n = nan(length(x),length(Name));
        for k = 1:length(Name)
            n(:,k) = IndexVal(Name{k},x);
        end
    end
    return
end
%
if isempty(Index)
    load('IndexData.mat','Index')   % Index structure
end
TabName = {Index.Name}; 
P = ismember(TabName,Name);
%
if nargin == 1
    n = Index(P).Fct;
else
    n = interp1(Index(P).n(:,1),Index(P).n(:,2),x);
    if isreal(Index(P).n(:,2)) && ~isempty(Index(P).k)
        n = n + 1i*interp1(Index(P).k(:,1),Index(P).k(:,2),x);
    end
end
%

end