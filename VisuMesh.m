function VisuMesh(Mesh,NumSD,str,alpha)

% VisuMesh
%   Plot 2D or 3D Mesh
%
% Syntax
%   VisuMesh(Mesh)
%   VisuMesh(Mesh,NumSD)
%   VisuMesh(Mesh,NumSD,str)
%   VisuMesh(Mesh,NumSD,str,alpha)
%
% Description
%   Mesh : Mesh of the structure (CoorN, Cn, Nsd,...CoorA, Ca, ExtAr, )
%   NumSD : Number of subdomais to be plotted,
%           if NumSD is a Matlab structure (Data ou Phys) we plot abs(Indice)
%   str   : edge color or 'none'
%   alpha : Transparency of colors for a 3D-mesh
%
%
% Example
%   h = [0.045 0.015 .05]; [dx,dy] = deal(.1,.1);
%   lix = {[] , 0.5*dx, 0};  liy = {[] , .75*dy, []};
%   npx = {2 [2 2] 2}; npy = {2 2 2}; npz = {2 2 2};
%   Mesh = MeshLayer(dx,lix,dy,liy,h,npx,npy,npz);
%   figure, VisuMesh(Mesh)
%   figure, VisuMesh(Mesh(2)), VisuMesh(Mesh([1 3]),1,'none',.3)

%
% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if nargin < 3 , str = 'k'; end
if nargin < 4 , alpha = 1; end

switch size(Mesh(1).CoorN,2)
    case 2
        if nargin == 1
            VisuMesh2D(Mesh);
        else
             VisuMesh2D(Mesh,NumSD,str);
        end
   case 3
        if nargin == 1
            VisuMesh3D(Mesh);
        else
             VisuMesh3D(Mesh,NumSD,str,alpha);
        end
end

end

%%
function VisuMesh2D(Mesh0,NumSD0,str)

%
% VisuMesh2D
%   Visualiser un Maillage 2D
%
% Syntaxe
%   VisuMesh2D(Mesh)
%   VisuMesh2D(Mesh,NumSD)
%   VisuMesh2D(Mesh,NumSD,str)
%   VisuMesh2D(Mesh,[],str)
%
% Description
%   Mesh : Stucture de données maillage
%   NumSD: Numéros des sous-domaines à afficher
%   str  : Couleur des arêtes ('k' par défaut, 'none' pas d'arêtes)
%
% Exemple
%   a1 = .8; l1 = [.1 .2 .1]; h = .2; npx = 8; npy = 10; ElmtOrder = 1;
%   Mesh = Mesh1Layer(a1,l1,h,npx,npy,ElmtOrder);
%   VisuMesh2D(Mesh)
%   VisuMesh2D(Mesh,[],'none')

% Date de la dernière version : 15 décembre 2017
% Auteur : Mondher Besbes (LCF / CNRS)

%figure
hold on
if nargin < 3 , str = 'k'; end

for kc = 1:length(Mesh0)
    Mesh = Mesh0(kc);

if nargin < 2
    [NumSD,C] = deal(1:max(Mesh.Nsd));
else
    if isempty(NumSD0)
        [NumSD,C] = deal(1:max(Mesh.Nsd));
    else
        if isstruct(NumSD0)
            Phys = NumSD0(kc);
            C = abs(Phys.Indice);
            NumSD = 1:length(Phys.Indice);
            str = 'none';
        else
            [NumSD,C] = deal(NumSD0);
        end
    end
end
%


%
for kd = 1:length(NumSD)
    k = NumSD(kd);
    Pe = find(Mesh.Nsd == k);
%
    [x,y] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2));
    [X,Y] = deal(x(Mesh.Cn(Pe,:)),y(Mesh.Cn(Pe,:)));
%
    fv = patch(X',Y',C(kd));
    set(fv,'FaceAlpha',1,'EdgeColor',str)
    axis equal,
    axis on, axis tight
    xlabel('x'), ylabel('y'),
    colorbar
end

% if length(Mesh0)>1,
%     plot([min(Mesh.CoorN(:,1)) max(Mesh.CoorN(:,1))],...
%         [max(Mesh.CoorN(:,2)) max(Mesh.CoorN(:,2))],'--r','LineWidth',4)
% end

end

end

%%
function VisuMesh3D(Mesh0,NumSD0,str,alpha)

if nargin < 2 , NumSD0 = []; end
if nargin < 3 , str = 'k'; end
if nargin < 4 , alpha = 1; end

%figure
if iscell(NumSD0)
    if length(Mesh0) == length(NumSD0)
        NumSD0 = SetData('Indice',NumSD0);
    elseif length(Mesh0) == length(NumSD0)-2
        NumSD0 = SetData('Indice',NumSD0(end-1:-1:2));
    end
end

NumMax = 0;

hold on
for kc = 1:length(Mesh0)
    Mesh = Mesh0(kc);
    if nargin < 2 || isempty(NumSD0),
        [NumSD,C] = deal(1:max(Mesh.Nsd));
    else
            if isstruct(NumSD0)
                Data = NumSD0(kc);
                Data = InterpIndex(Data);
                C = abs(Data.Indice);
                NumSD = 1:length(Data.Indice);
            else
                [NumSD,C] = deal(NumSD0);
            end
    end
    if kc>=2 && isfield(Mesh,'dh'), Mesh = MoveMesh(utilMesh(Mesh),[0 0 Mesh.dh]); end

    TabPf = 1:size(Mesh.ExtFa,1);

    for ki = 1:length(NumSD),
        k = NumSD(ki);

        Pe = Mesh.Nsd == k;
        if isfield(Mesh,'Pe') && ~isempty(Mesh.Pe), Pe = Pe & eval(char(Mesh.Pe)); end
        Pf = unique(Mesh.Cf(Pe,:));
        if isfield(Mesh,'CoorA') && k<=max(Mesh.Nsd)
            if isfield(Mesh,'TabNsdF'), Pf = find(ismember(TabPf(:),Pf) & Mesh.TabNsdF{k} == k); end
        else
            Pf = find(ismember(TabPf(:),Pf) & Mesh.TabP == 1);
        end
%        
        if alpha < 1, Pf = find(Mesh.TabNsdF{k} == k); end %Pf = find(Mesh.TabP == 1); end

%        Label{k} = num2str(C(ki));
        Label{kc} = num2str(max(C));

        [x,y,z] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2),Mesh.CoorN(:,3));
        [X,Y,Z] = deal(x(Mesh.ExtFa(Pf,:)),y(Mesh.ExtFa(Pf,:)),z(Mesh.ExtFa(Pf,:)));


        fv = patch(X',Y',Z',C(ki));
        set(fv,'FaceAlpha',alpha,'EdgeColor',str)
        axis equal,
        axis on, view(3), axis tight
        xlabel('x'), ylabel('y'), ylabel('y')
        drawnow
        if length(Mesh0) ==1 && nargin >= 2 && k < max(NumSD), pause(1), end

    end
%
if ~(isfield(Mesh,'Pe') && ~isempty(Mesh.Pe)), 
if isfield(Mesh,'xv') && ~isempty(Mesh.xv)
    for k = 1:size(Mesh.xv,1)
    plot3(Mesh.xv(k,:),Mesh.yv(k,:),max(Mesh.CoorN(:,3))*ones(size(Mesh.xv(k,:))),'k','LineWidth',1);
    plot3(Mesh.xv(k,:),Mesh.yv(k,:),min(Mesh.CoorN(:,3))*ones(size(Mesh.xv(k,:))),'k','LineWidth',1);
    end
end
end

NumMax = ceil(max([NumMax Label{kc}]));
%colorbar('Ticks',NumSD,'TickLabels',Label)
box on,
set(gca,'BoxStyle','full','LineWidth',2);

[xm,xM,ym,yM] = MinMax(util(Mesh.CoorN(:,1:2)));
axis([xm,xM,ym,yM])

end
%colorbar('Ticks',NumSD,'TickLabels',Label)
colorbar('Ticks',1:NumMax)

end






