function Data = SetData(varargin)

% 
% SetData 
%   Define data of light sources,  physical properties, ...
%
% Syntax
%   Data = SetData('Param1',Val1,'Param2',Val2,...);
%   Data = SetData(); % Default value
%
% Description
%   Param     : Name of parameter
%   Val       : Value of Param
%   Geom.Help : Parameter lexicon of the structure Data

% Date of the latest version : 06 february 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

          
%if nargin == 0, return, end
if nargin == 1, Data = varargin{1}; return, end
DataIn = 0;
if nargin > 1 & isstruct(varargin{1})
    DataIn = 1;
    Data = varargin{1}; 
    varargin = varargin(2:end); 
else
    Data = InitSetData;
end

% Mise à jour des champs de la structure Data
% -------------------------------------------

f = varargin(1:2:end); % nom des champs à modifier
c = varargin(2:2:end); % valeurs associées
%
for kc = 1:length(Data)
    FieldData = fieldnames(Data(kc));

    for k = 1:length(f)
        P = ismember(FieldData,f(k));
        if any(P) 
            Data(kc).(FieldData{P}) = c{k};
        else % nouveau champ
            Data(kc).(f{k}) = c{k};
            warning(['Nouvelle donnée : ', f{k}]); 
        end
    end
end

if DataIn == 0
    P1 = find(ismember(f,'Indice'));
    P2 = find(ismember(f,'TabIndice'));
    if isempty(P2) 
        Data.TabIndice = Data.Indice; 
    elseif isempty(P1) 
        Data.Indice = InterpIndex(Data.TabIndice,Data.Lambda0); 
    end
end


Indice1 = Data.Indice;
%
if iscell(Indice1)  %&& length(Data) == 1
    Data0 = Data;
    clear Data
    for kc = 1:length(Indice1)
%        if isfield(Data0(1),'my')  && ~isempty(Data0(1).my), 
            Data(kc) = SetData(Data0(1),'Indice',nan(1,length(Indice1{kc})));
            Data(kc).TabIndice = Indice1{kc};
            if iscell(Data(kc).TabIndice)
                for k = 1:length(Data(kc).TabIndice)
                    Data(kc) = InterpIndex(Data(kc));
                end
            else
                Data(kc).Indice = Indice1{kc};
            end

    end
    return
end
% %

%
for kc = 1:length(Data)

if isfield(Data(kc),'my') && ~isempty(Data(kc).my), return, end
if isfield(Data(kc),'my') && isempty(Data(kc).my) 
    if  Data(kc).TypePol == 0
        Data(kc).DofType = 'Node';
    end
end


end
end



%%
function Data = InitSetData()

% Lexique des champs
% ------------------
Lexique = struct('Lambda0','Wavelength (µm)',...
              'Theta0', 'Angle of incidence (rd)',...
              'Phi0','Azimuthal angle  (rd)',...
              'Psi0','Polarization angle  (rd)',...
              'a','Xmax to define PML position in x, by default [-Inf +Inf]',...
              'b','Ymax to define PML position in y, by default [-Inf +Inf]',...
              'c','Zmax to define PML position in z, by default [-Inf +Inf]',...
              'Rmax','Rmax to define Dirichlet condition in spherical surface',...
              'Indice','Refractive index of different medium',...
              'nb','Refractice index of substrat',...
              'nh','Refractice index of supstrat',...
              'TypePol','0: TE  - 2: TM, by default =2',...
              'CoefPml','Coefficient of PML',...
              'Ie','Electric Dipole, by default [0;0;0]',...
              'Im','Magnetic Dipole, by default [0;0;0]',...
              'CoorDip','Coordonnate of the dipole',...
              'ChampInc',' = -1: from down, +1: from top',... 
              'Beam','LGB ; GB ; PW ; Dipole',...
              'x0','Position of the waist in x',...
              'y0','Position of the waist in y',...
              'z0','Position of the waist in z',...
              'w0','Waist of the beam (µm)',...
              'mx','Number of Fourier terms in x',...
              'my','Number of Fourier terms y',...
              'Sym','0: Neumann  1: Dirichlet 2: Periodic',...
              'DofType','DOF Node ou Edge');
%
Lexique_Fr = struct('Lambda0','Longueur d''onde en µm',...
              'Theta0', 'Angle d''incidence',...
              'Phi0','Angle d''incidence',...
              'Psi0','Angle d''incidence',...
              'a','Xmax pour définir les PML',...
              'b','Ymax',...
              'c','Zmax',...
              'Rmax','Rmax',...
              'Indice','Indice des différents milieux',...
              'nb','Indice du bas',...
              'nh','Indice du haut',...
              'TypePol','0: TE  - 2: TM',...
              'CoefPml','Coefficient du PML',...
              'Ie','Dipôle électrique',...
              'Im','Dipôle magnétique',...
              'CoorDip','Coordonnées du dipôle',...
              'ChampInc','Champ incident',... 
              'Beam','LGB ; GB ; PW ; Dipole',...
              'x0','Position du faisceau',...
              'y0','Position du faisceau',...
              'z0','Position du faisceau',...
              'w0','west',...
              'mx','Nombre de termes de Fourier en x',...
              'my','Nombre de termes de Fourier en y',...
              'Sym','0: Neuman  1: Dirichlet 2: Périodicité',...
              'Period','Période du réseau',...
              'hc','Hauteur des couches',...
              'nc','Indice des couches avec FMM',...
              'li','Largeur des inclusions',...
              'DofType','DOF Node ou Edge');

% Valeurs par défaut
% ------------------
Data = struct('Lambda0', .6 , ...       % Longueur d'onde en µm
              'Theta0', 0, ...          % Angle d'incidence
              'Phi0', 0, ...            % Angle d'incidence
              'Psi0', 0, ...            % Angle d'incidence
              'a',[-Inf +Inf], ...      % Xmax pour définir les PML
              'b',[-Inf +Inf], ...      % Ymax
              'c',[-Inf +Inf], ...      % Zmax
              'Rmax',Inf, ...           % Rmax
              'Indice',1,...            % Indice des différents milieux
              'nb',1,...                % Indice du bas
              'nh',1,...                % Indice du haut
              'TypePol',2, ...          % 0: TE  - 2: TM
              'CoefPml',5*(1+1i), ...   % Coefficient du PML
              'Ie',[0; 0; 0], ...       % Dipôle électrique
              'Im',[0; 0; 0], ...       % Dipôle magnétique
              'CoorDip',[0 0 0],...     % Coordonnées du dipôle
              'ChampInc',+1,...         % Champ incident 
              'Beam','PW',...           % LGB ; GB ; PW ; Dipole
              'x0',0.00,...             % Position du faisceau
              'y0',0.00,...             % Position du faisceau
              'z0',0.00,...             % Position du faisceau
              'w0',1,...                % west
              'mx',0,...                % Nombre de termes de Fourier en x
              'my',[],...                % Nombre de termes de Fourier en y
              'Sym',[2 2],...           % 0: Neuman  1: Dirichlet 2: Périodicité
              'Period',.2,...           % Période du réseau
              'hc',0,...               % Hauteur d'une couche
              'nc',1,...                % Indice d'une couche
              'li',0,...                % Largeur des inclusions
              'DofType','Edge',...      % DOF Node ou Edge
              'Help',Lexique); 

end
          
