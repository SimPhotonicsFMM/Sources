function [Geom,Mesh] = SetGeom(varargin)

% 
% SetGeom
%   Define geometric parameters and plot the shape of particle defined by
%   superformula if parameter 'Plot' is equal to 1.
%
% Syntax
%   Geom = SetGeom('Param1',Val1,'Param2',Val2,...);
%   Geom = SetGeom(); % Default value
%
% Description
%   Param     : Name of parameter
%   Val       : Value of Param
%   Geom.Help : Parameter lexicon of the structure Geom
%
% See also TutoGeomMesh.mlx 

% Date of the latest version : 21 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)


% Lexique des champs
% ------------------
Lexique = struct('hc','Layer thickness [hc_1 ..hc_i..] (µm)',...
        'dx','Period in x-axis (µm)',...
        'dy','Period in y-axis (µm)',...
        'mn','Superformula parameters [m, n1, n2, n3]',...
        'ab','Radius [a b] respectively in x and y (µm)',...
        'Angle','Angle of inclusion rotation (rd)',...
        'Dep','Translation in x and y (µm)',...
        'NumSD','Subdomain number of inclusions, otherwise =1 for all inclusions',...
        'Np','Number of points for inclusion drawing (Np = 401 by default)',...
        'npx','Number of nodes in x (by default npx=2 as min value per inclusion)',...
        'npy','Number of nodes in y (as npx)',...
        'npz','Number of nodes in z, by default npz=2',...
        'Plot','=1: Plot particules, 0: no plot');
        
%        'lix','Width Inclusions in x, by default lix=[]',... 
%        'liy','Width Inclusions in y, by default liy=[]',...
%        'R','Radius of inclusions',...

%
% French lexicon
Lexique_Fr = struct('dx','Période en x (µm)',...
        'dy','Période en y (µm)',...
        'hc','Tableau hauteurs des différentes couches [hc_1 .. hc_i...] en µm',...
        'lix','Largeur Inclusions en x, par défaut lix=[]',... 
        'liy','Largeur Inclusions en y, par défaut liy=[]',...
        'npx','Nombre de noeuds en x (=2 ou valeur min par inclusion)',...
        'npy','Nombre de noeuds en y',...
        'npz','Nombre de noeuds en z, par défaut npz=2',...
        'R','Rayons des inclusions',...
        'mn','Liste des paramètres [m, n1, n2, n3] de la superformule',...
        'ab','Rayons [a b] respectivement en x et y en µm',...
        'Angle','Angle de rotation en rd',...
        'Dep','Translation en x et y en µm',...
        'Np','Nombre de points (Np = 401 par défaut)',...
        'NumSD','Numéro des sous-domaines des nano-particules',...
        'Plot','=1:Tracer les contours des particules, 0:sinon',...
        'Echelle','Echelle de la géométrie, par défaut =1 (en µm)', ...
      	'TabSym','Tableau des symétries 1: x=0 , 2: y=0 , 3: z=0 , []: pas de symétrie',...
        'Sym','Changer (1) ou non (0) les numéros des sous-domaines après une symétrie',...
    	'ElmtOrder','Ordre des éléments',...
        'Dim','Dimension du problème 2 ou 3',...
        'File','Nom du fichier Géométrie pour le mailleur Cubit');
    
% Valeurs par défaut
% ------------------
Geom = struct('Help',Lexique, ...
              'ElmtOrder',1,...
              'TabSym',[],...
              'Sym',0,...
              'Echelle',1,...
              'Dim',3);
          
if nargin == 1, Geom = varargin{1}; return, end
if nargin > 1 && isstruct(varargin{1}) 
    Geom = varargin{1}; 
    varargin = varargin(2:end); 
end

% Mise à jour des champs de la structure Data
% -------------------------------------------

f = varargin(1:2:end); % nom des champs à modifier
c = varargin(2:2:end); % valeurs associées
%
FieldData = fieldnames(Geom);

for k = 1:length(f)
    P = ismember(FieldData,f(k));
    if any(P) 
        Geom.(FieldData{P}) = c{k};
    else % nouveau champ
        Geom.(f{k}) = c{k};
        %warning(['Nouvelle donnée : ', f{k}]); 
    end
end
%
if nargin == 0, return; end
% 0D grating
if ~isfield(Geom,'dx') && ~isfield(Geom,'dy'), [Geom.dx,Geom.dy] = deal(max(Geom.hc)); end
% 1D grating
if isfield(Geom,'dx') && ~isfield(Geom,'dy'), Geom.dy = Geom.dx; end
%
if isfield(Geom,'ab')
    if iscell(Geom.ab)
        for kc = 1:length(Geom.ab)
            if isempty(Geom.ab{kc}), Geom.ab{kc}=[Geom.dx/2 Geom.dy/2]; end
            P = isinf(Geom.ab{kc}(:,2)) | isnan(Geom.ab{kc}(:,2));
            Geom.ab{kc}(P,2) = Geom.dy/2;
        end
    else
        if isempty(Geom.ab), Geom.ab=[Geom.dx/2 Geom.dy/2]; end
        P = isinf(Geom.ab(:,2)) | isnan(Geom.ab(:,2));
        Geom.ab(P,2) = Geom.dy/2;
    end
end
%
if nargout == 2
    Mesh = MeshLayer(Geom);
end
% 
if isfield(Geom,'Plot') && Geom.Plot == 1
    %
    c = {'b' 'g' 'r' 'k' 'y' 'c'};
    if isfield(Geom,'mn')
        % 
        figure, hold on
        %
        [xv,yv] = deal([]);
        if iscell(Geom.mn)
            for kc = 1:length(Geom.mn)
                n0 = Geom.mn{kc};
                if size(n0,2) == 1, n0 = [n0 repmat([400 400 400],size(n0,1),1)]; end
                a0 = Geom.ab{kc};
                if isfield(Geom,'Np'), Np = Geom.Np; else, Np = 401; end
                if isfield(Geom,'Dep'), Dep0 = Geom.Dep{kc}; else, Dep0 = zeros(size(a0)); end
                if isfield(Geom,'Angle'), angle0 = Geom.Angle{kc}; else, angle0 = zeros(size(n0,1),1); end
    
                if iscell(Np), Np = cell2mat(Np); end
                [x0,y0] = sf2d(n0, a0, angle0, Dep0, Np);
                xv = [xv ; x0];
                yv = [yv ; y0];
            end
            
        else
            [xv,yv] = sf2d(Geom);
        end
        %
        for k = 1:size(xv,1)
            PlotCourbe(xv(k,:),yv(k,:),'','','',c{randi(length(c))})
        end
        %
        axis equal
        axis([-Geom.dx/2 Geom.dx/2 -Geom.dy/2 Geom.dy/2])
    end
end
%


return


          