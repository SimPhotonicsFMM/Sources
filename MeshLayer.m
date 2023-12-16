function Mesh = MeshLayer(varargin)

global SaveOption

% MeshLayer 
%   Generate a 2D mesh in quadrangles or 3D in parallelelipeds
%       SaveOption = 1 all fields of the structure "Mesh" otherwise = 0
%
% Syntax
%   Mesh = MeshLayer(Geom);                   
%   Mesh = MeshLayer(Geom,CoefFlip);                   
%   Mesh = MeshLayer(a1,l1,h,npx,npy);                  % 2D
%   Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz);        % 3D
%   Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz,R); 
%   Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz,xv,yv);   
%   Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz,xv,yv,NumSD);   
%
% Description
%   Geom : Structure genearated by "SetGeom' of geometric parameters 
%   CoefFlip : =1 to flip layers in up/down direction
%   a1 : Period of the structure (on the x-axis)
%   a2 : Period of the structure (on the y-axis)
%   h  : Layer thickness
%
% Array or cell array if multiple layers
%   l1 : Width Inclusions in x, l1=[] no inclusion
%   l2 : Width Inclusions in y, l2=[] no inclusion
%   npx : Number of nodes in x (or npx=2 as min value per inclusion)
%   npy : Number of nodes in y
%   npz : Number of nodes in z
%   R   : Radius of inclusions
%   xv  : x-coordinates of points describing an inclusion
%   yv  : y-coordinates of points describing an inclusion
%   NumSD : Subdomain number of inclusions
%
%   Mesh : Stucture array for each layer
%       Cn	: Connectivity to nodes
%       Ca	: Connectivity to edges
%       Cf	: Connectivity to facets
%       CoorN	: Node coordinates
%       CoorA	: Edge coordinates
%       CoorF	: Facte coordinates
%       CoorV	: Volume coordinates
%       Nsd	: Subdomain numbers
%       ExtAr	: Numbers of the nodes ends of the edges
%       ExtFa	: Numbers of the nodes ends of the facets 
%       ExtFaAr	: Numbers of the edges ends of the facets 
%       TypeElmt: 'Hexahedron'
%       NsdF	: Numbers of subdomains associated with the facets
%       TabP	: Pointer to outer facets 
%       Surf	: Surface of facets 
%       Vol	: Volume of mesh elements
%       TabNsdF: Numbers of subdomains associated with the facets per subdomain
%
% Example 2D
%   a1 = 1; l1 = [.1 .2 .3]; h = .2; npx = 8; npy = 10; 
%   Mesh = MeshLayer(a1,l1,h,npx,npy);
%
% Example 3D
%   h = [0.045 0.015]; [dx,dy] = deal(.1,.1);      
%   lix = {[] , 0.5*dx};  liy = {[] , .75*dy};
%   npx = {2 [2 2]}; npy = {2 2}; npz = {2 2};
%   Mesh = MeshLayer(dx,lix,dy,liy,h,npx,npy,npz);
%   figure, VisuMesh(Mesh)
%
% See also TutoGeomMesh.mlx, SetGeom.m

% Date of the latest version : 24 January 2022
% Author : Mondher Besbes (LCF / CNRS / IOGS)
% 

if isstruct(varargin{1,1})
    Geom = varargin{1,1};
    Geom = FlipUD(util(Geom));
    if isfield(Geom,'SaveOption'), SaveOption = Geom.SaveOption; else SaveOption = 1; end
    %
    if length(varargin) == 2 && varargin{1,2} == 1, Geom = FlipUD(util(Geom)); end
    %    
    if length(varargin) >= 3
        [npx,npy] = deal(varargin{1,2},varargin{1,3});
    else
        if isfield(Geom,'npx')
            [npx,npy] = deal(Geom.npx,Geom.npy);
        else
            [npx,npy] = deal(2);
        end
    end
    %
    [a1,h] = deal(Geom.dx,Geom.hc);
    if ~isfield(Geom,'dy') % cas geométrie 2D
        if isfield(Geom,'lix')
            l1 = Geom.lix;
        else
            if length(h)>1, l1 = cell(length(h)); else, l1 = []; end
        end

        Mesh = MeshLayer(a1,l1,h,npx,npy);
    else
        if length(varargin) == 4
            npz = deal(varargin{1,4});
        else
            if isfield(Geom,'npz')
                npz = Geom.npz;
            else
                npz = 2;
            end
        end
        a2 = Geom.dy;
        if isfield(Geom,'lix') && isfield(Geom,'liy')
            [l1,l2] = deal(Geom.lix,Geom.liy);
        else
            if length(h)>1, [l1,l2] = deal(cell(length(h))); else, [l1,l2] = deal([]); end
        end
            
        if isfield(Geom,'R')
            Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz,Geom.R);
        elseif isfield(Geom,'mn') 
            if iscell(Geom.mn)
                [l1,l2] = deal(cell(1,length(h)));
                for k = 1:length(Geom.mn)
                    Geom.Plot = 0;
                    Geom1 = SetGeom(Geom,'mn',Geom.mn{k},'ab',Geom.ab{k});
                    if isfield(Geom,'Np'), Geom1.Np = Geom.Np{k}; end
                    if isfield(Geom,'Angle')
                        Geom1.Angle = Geom.Angle{k};
                    else
                        Geom1.Angle = zeros(1,size(Geom1.mn,1));
                    end
                    if isfield(Geom,'Dep')
                        Geom1.Dep = Geom.Dep{k};
                    else
                        Geom1.Dep = zeros(size(Geom1.mn,1),2);
                    end
                    if isfield(Geom,'NumSD'), NumSD{k} = Geom.NumSD{k}; else, NumSD = []; end
                    Geom1.Plot = 0;
                    [xv{k},yv{k}] = sf2d(Geom1);
                    if isempty(l1{k}) && isempty(l2{k})
                        if (size(Geom1.mn,2) == 4 || any(Geom1.Angle ~= 0))...
                                || (max2(util(npx))>2 || max2(util(npy))>2)
                            l1{k} = max(xv{k}(:))-min(xv{k}(:));
                            l2{k} = max(yv{k}(:))-min(yv{k}(:));
                            %
                            if l1{k}>=Geom1.dx, l1{k} = []; end
                            if l2{k}>=Geom1.dy, l2{k} = []; end
                        else
                            l1{k} = diff(sort(unique2(util(xv{k}(:)),1e-10)));
                            l2{k} = diff(sort(unique2(util(yv{k}(:)),1e-10)));
                            l1{k}(l1{k}<10*eps)=[];
                            l2{k}(l2{k}<10*eps)=[];

                        end
                    end
                end
            else
                if isfield(Geom,'NumSD'), NumSD = Geom.NumSD; else, NumSD = []; end
                if ~isfield(Geom,'Angle')
                   Geom.Angle = zeros(1,size(Geom.mn,1));
                end
                Geom.Plot = 0;
                [xv,yv] = sf2d(Geom);
                if isempty(l1) && isempty(l2)
                    if size(Geom.mn,2) == 4 || any(Geom.Angle ~= 0) || ...
                            max(npx)>2 ||  max(npy)>2
                        l1 = max(xv(:))-min(xv(:));
                        l2 = max(yv(:))-min(yv(:));
                        %
                        if l1>=Geom.dx, l1 = []; end
                        if l2>=Geom.dy, l2 = []; end
                    else
                        l1 = diff(sort(unique2(util(xv(:)),1e-10)));
                        l2 = diff(sort(unique2(util(yv(:)),1e-10)));
                        l1(l1<10*eps)=[];
                        l2(l2<10*eps)=[];
                    end
                end
            end
            Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz,xv,yv,NumSD);

        else
            Mesh = MeshLayer(a1,l1,a2,l2,h,npx,npy,npz);
        end
    end
    return
end
%
[a1,l1] = deal(varargin{1,1},varargin{1,2});
if size(varargin,2) == 5
    [h,npx,npy] = deal(varargin{1,3},varargin{1,4},varargin{1,5});
    Mesh = MeshLayer2D(a1,l1,h,npx,npy);
else
    [a2,l2] = deal(varargin{1,3},varargin{1,4});
    [h,npx,npy,npz] = deal(varargin{1,5},varargin{1,6},varargin{1,7},varargin{1,8});
    if length(h)>1 || iscell(l1)
        dh=cumsum(h);
        dh = [0 dh(1:end-1)];
        for kc = 1:length(h)
            if iscell(l1(kc)), l1c = l1{kc}; else, l1c = l1; end
            if iscell(l2(kc)), l2c = l2{kc}; else, l2c = l2; end
            if iscell(npx), npxc = npx{kc}; else, npxc = npx; end
            if iscell(npy), npyc = npy{kc}; else, npyc = npy; end
            if iscell(npz), npzc = npz{kc}; else, npzc = npz; end
            %
            if size(varargin,2) == 8
                Mesh(kc) = MeshLayer3D(a1,l1c,a2,l2c,h(kc),npxc,npyc,npzc);
            elseif size(varargin,2) == 9
                R = varargin{1,9};
                if ~isempty(R) && iscell(R(kc)), Rc = R{kc}; else, Rc = R; end
                if Rc == 0, Rc = []; end
                Mesh(kc) = MeshLayer3D(a1,l1c,a2,l2c,h(kc),npxc,npyc,npzc,Rc);
            else
                xv = varargin{1,9};
                yv = varargin{1,10};
                if iscell(xv), xvc = xv{kc}; else, xvc = xv; end
                if iscell(yv), yvc = yv{kc}; else, yvc = yv; end
                if size(varargin,2) == 11
                    NumSD = varargin{1,11};
                    if iscell(NumSD), NumSDc = NumSD{kc}; else, NumSDc = NumSD; end
                    Mesh(kc) = MeshLayer3D(a1,l1c,a2,l2c,h(kc),npxc,npyc,npzc,xvc,yvc,NumSDc);
                else
                    Mesh(kc) = MeshLayer3D(a1,l1c,a2,l2c,h(kc),npxc,npyc,npzc,xvc,yvc);
                end
            end
            if SaveOption == 1
                Mesh(kc) = MoveMesh(utilMesh(Mesh(kc)),[0 0 dh(kc)]);
            else
                Mesh(kc) = MoveMesh(utilMesh2(Mesh(kc)),[0 0 dh(kc)]);
            end
        end
    else

        if size(varargin,2) == 8
            Mesh = MeshLayer3D(a1,l1,a2,l2,h,npx,npy,npz);
        elseif size(varargin,2) == 9
            R = varargin{1,9};
            Mesh = MeshLayer3D(a1,l1,a2,l2,h,npx,npy,npz,R);
        else
            xv = varargin{1,9};
            yv = varargin{1,10};
            if size(varargin,2) == 11
                NumSD = varargin{1,11};
                Mesh = MeshLayer3D(a1,l1,a2,l2,h,npx,npy,npz,xv,yv,NumSD);
            else
                Mesh = MeshLayer3D(a1,l1,a2,l2,h,npx,npy,npz,xv,yv);
            end

        end
    end
end
    
end




%%
function Mesh = MeshLayer2D(a1,l1,h,npx,npy)

% MeshLayer2D 
%   Generate a 2D mesh in quadrangles

%
if nargin < 4, npx = 2; npy = 2; end
%
if length(h)>1
    dh=cumsum(h);
    dh = [0 dh(1:end-1)];
    for kc = 1:length(h)
        if iscell(l1(kc)), l1c = l1{kc}; else, l1c = l1(kc); end
        if iscell(npx), npxc = npx{kc}; else, npxc = npx; end
        if iscell(npy), npyc = npy{kc}; else, npyc = npy; end
        Mesh(kc) = MeshLayer(a1,l1c,h(kc),npxc,npyc);
        Mesh(kc) = MoveMesh(utilMesh(Mesh(kc)),[0 dh(kc)]);
    end
    return
end

% Valeur npx si plusieurs inclusions
if l1 == 0, l1 = []; end
if length(npx) == 1
if npx > 3
    dx = min([l1(:);a1-sum(l1)])/(npx);
    npx = ceil([a1/2-sum(l1)/2; l1(:); a1/2-sum(l1)/2]/dx);
else
    npx = repmat(npx,length(l1)+2,1);
end
end

% Coordonnées en x et y
if isempty(l1) || all(l1 == 0)
    x0 = linspace(0,a1,npx(1));
else
    x0 = linspace2(util(0),util(a1/2-sum(l1)/2),npx(1),1,1);
    xi = a1/2-sum(l1)/2;
    for k = 1:length(l1)
        x0 = [x0, linspace2(util(xi),util(xi+l1(k)),npx(k+1),1,1)];
        xi = xi+l1(k);
    end
    x0 = [x0, linspace2(util(xi),util(a1),npx(end),1)]; 
end
x0 = x0-a1/2;
%
y0 = linspace(0,h,npy);

% Maillage

[x,y] = meshgrid(x0,y0);

Cn = [];
for i = 1:length(y0)-1
    i0 = (1:length(x0)-1) + (i-1)*length(x0);
    i1 = i0 + length(x0);
    Cn = [Cn; [i0' i0'+1 i1'+1 i1']];
end
    
CoorN = [];
j = 1:length(x0);
for i = 1:length(y0)
    CoorN = [CoorN; [x(i,j)' y(i,j)']];
end

% Numéro des sous-domaines

Nsd = (length(l1)+1)*ones(size(Cn,1),1);    % domaine extrême

[x,y] = deal(CoorN(:,1),CoorN(:,2));
[X,Y] = deal(x(Cn),y(Cn));
[xg,yg] = deal(sum(X,2)/4,sum(Y,2)/4);

xi = -sum(l1)/2;
for k = 1:length(l1)
    Nsd(xg>xi & xg<xi+l1(k)) = k; 
    xi = xi+l1(k);
end

% Recherche des numéros des arêtes ...
[Ca,ExtAr,CoorA,CoorF] = RechercheArete(Cn,CoorN);
TypeElmt = 'Quadrangle';

% Structure Mesh
Mesh = struct('Cn',Cn,'Ca',Ca,'CoorN',CoorN,'CoorA',CoorA,'CoorF',CoorF,...
              'Nsd',Nsd,'ExtAr',ExtAr,'TypeElmt',TypeElmt);

end

%%
function Mesh = MeshLayer3D(a1,l1,a2,l2,h,npx,npy,npz,R,yv,NumSD)

global SaveOption

% MeshLayer3D 
%   Generate a 3D mesh in parallelelipeds

% 
%
if nargin < 11, NumSD = []; end
if nargin < 10, xv = []; yv = []; end
if nargin < 9, R = []; end
if nargin < 5, npx = 2; npy = 2; npz = 2; end

if nargin > 9, xv = R; R = []; end

%

% Valeur npx si plusieurs inclusions
if l1 == 0, l1 = []; end
if l2 == 0, l2 = []; end
%
if length(npx) ==1
    npx = repmat(npx,length(l1)+2,1); 
else
    npx = [npx(end) npx];
end
%
if length(npy) ==1 
    npy = repmat(npy,length(l2)+2,1); 
else
    npy = [npy(end) npy];
end

% Coordonnées en x et y
if isempty(l1) || all(l1 == 0)
    x0 = linspace(-a1/2,a1/2,npx(1));
else
    if ~isempty(xv)
        x0 = linspace2(util(-a1/2),util(min(xv(:))),npx(1),1,1);
        if abs(-a1/2-min(xv(:)))<10*eps, x0 = []; end
        xi = min(xv(:));
        if -a1/2>min(xv(:)), x0 = []; xi = -a1/2; end
    else
        x0 = linspace2(util(-a1/2),util(-sum(l1)/2),npx(1),1,1);
        if abs(-a1/2+sum(l1)/2)<10*eps, x0 = []; end
        xi = -sum(l1)/2;
    end
    for k = 1:length(l1)
        x0 = [x0, linspace2(util(xi),util(xi+l1(k)),npx(k+1),1,1)];
        xi = xi+l1(k);
        if abs(xi)<10*eps, xi = 0; end
    end
    if abs(a1/2-xi)<10*eps
        x0 = [x0,a1/2];
    else
        x0 = [x0, linspace2(util(xi),util(a1/2),npx(end),1)];
    end
end
x0 = unique2(util(sort([x0, min(xv,[],2)', max(xv,[],2)'])),1e-10);
x0 = x0(abs(x0)<=a1/2+eps);

%x0 = x0-a1/2;
%
if isempty(l2) || all(l2 == 0)
    y0 = linspace(-a2/2,a2/2,npy(1));
else
    if ~isempty(yv)
        y0 = linspace2(util(-a2/2),util(min(yv(:))),npy(1),1,1);
        if abs(-a2/2-min(yv(:)))<10*eps, y0 = []; end
        yi = min(yv(:));
        if -a2/2>min(yv(:)), y0 = []; yi = -a2/2; end
    else
        y0 = linspace2(util(-a2/2),util(-sum(l2)/2),npy(1),1,1);
        if abs(-a2/2+sum(l2)/2)<10*eps, y0 = []; end
        yi = -sum(l2)/2;
    end
    for k = 1:length(l2)
        y0 = [y0, linspace2(util(yi),util(yi+l2(k)),npy(k+1),1,1)];
        yi = yi+l2(k);
        if abs(yi)<10*eps, yi = 0; end
    end
    if abs(a2/2-yi)<10*eps
        y0 = [y0,a2/2];
    else
        y0 = [y0, linspace2(util(yi),util(a2/2),npy(end),1)];
    end
end
y0 = unique2(util(sort([y0, min(yv,[],2)', max(yv,[],2)'])),1e-10);
y0 = y0(abs(y0)<=a2/2+eps);
%y0 = y0-a2/2;

if nargin <= 11
    if ~isempty(yv) && npz==2, [x0,y0,dx,dy] = OptimMesh(x0,y0,xv,yv,NumSD); end
end
%
z0 = linspace(0,h,npz);

% Maillage
P0 = find(abs(x0)<eps);
if isempty(P0), x0 = sort([0 x0]); end
P0 = find(abs(y0)<eps);
if isempty(P0), y0 = sort([0 y0]); end

    
[x,y] = meshgrid(x0,y0);

% Maillage de la surface de base
% ------------------------------

% Cn = [];
% for in = 1:length(y0)-1
%     i0 = (1:length(x0)-1) + (in-1)*length(x0);
%     i1 = i0 + length(x0);
%     Cn = [Cn; [i0' i0'+1 i1'+1 i1']];
% end
Cn = nan((length(x0)-1)*(length(y0)-1),4);%
%Cn = [];
for in = 1:length(y0)-1
    i0 = (1:length(x0)-1) + (in-1)*length(x0);
    i1 = i0 + length(x0);
    %Cn = [Cn; [i0' i0'+1 i1'+1 i1']];
    Cn((1:length(x0)-1)+(in-1)*(length(x0)-1),:) = [i0' i0'+1 i1'+1 i1'];
end
    
% CoorN = [];
% jn = 1:length(x0);
% for in = 1:length(y0)
%     CoorN = [CoorN; [x(in,jn)' y(in,jn)' z0(1)*ones(size(y(in,jn)))']];
% end

%CoorN = [];
CoorN = nan((length(x0))*(length(y0)),3);
jn = 1:length(x0);
for in = 1:length(y0)
    %CoorN = [CoorN; [x(i,j)' y(i,j)']];
    CoorN(jn+(in-1)*length(x0),:) = [x(in,jn)' y(in,jn)' z0(1)*ones(size(y(in,jn)))'];
end


% Maillage 3D
% -----------

CoorN0 = CoorN;
for k = 2:length(z0)
    CoorN = [CoorN; [CoorN0(:,1:2) CoorN0(:,3)+z0(k)-z0(1)]];
end

Cn0 = Cn;
ExtFa = Cn;
for k = 2:length(z0)
    ExtFa = [ExtFa; Cn0+(k-1)*length(CoorN0)];
end

Cn = [];
for k = 1:length(z0)-1
    Cn = [Cn; [ExtFa((k-1)*size(Cn0,1)+1:k*size(Cn0,1),:),ExtFa(k*size(Cn0,1)+1:(k+1)*size(Cn0,1),:)]];
end

% Recherche des arêtes
% --------------------

Na = 12;
% Initialisation
ExtAri = zeros(Na*size(Cn0,1),2);            % Extrémités arêtes
Cai = [];                                   % Connectivités arêtes
CoorV = zeros(size(Cn0,1),size(CoorN,2));    % Coordonnées des volumes

Tab = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
for k = 1:1:Na, ExtAri(k:Na:Na*size(Cn,1),:) = Cn(:,Tab(k,:)); end
for k = 1:1:Na, Cai = [Cai, Na*(1:size(Cn,1))'-(Na-k)]; end

% Eliminer les doublons
[ExtAr,I,J] = unique(sort(ExtAri,2),'rows');
ExtAr = ExtAri(I,:);
Ca = J(Cai);
if size(Ca,2) == 1, Ca = Ca'; end

% Coordonnées des arêtes
TabCoorA = zeros(length(ExtAr),size(CoorN,2),size(ExtAr,2)); 
for k = 1:1:size(TabCoorA,3)
    TabCoorA(:,:,k) = CoorN(ExtAr(:,k),:);
end
CoorA = sum(TabCoorA,3)/size(TabCoorA,3);

                  
% Signes des arêtes et coordonnées des volumes
TabCoorV = zeros(size(Cn,1),size(CoorN,2),size(Cn,2)); 
for k = 1:1:size(TabCoorV,3)
    TabCoorV(:,:,k) = CoorN(Cn(:,k),:);
end
CoorV = sum(TabCoorV,3)/size(TabCoorV,3);
%
for ie = 1:size(Ca,1),
    for ia = 1:1:size(Ca,2), if ExtAr(Ca(ie,ia),1) ~= Cn(ie,Tab(ia,1)), Ca(ie,ia) = -Ca(ie,ia); end, end
    %CoorV(ie,:) = sum(CoorN(Cn(ie,:),:),1)/size(Cn,2);
end

% Recherche des facettes
% ----------------------

%tic
Nf = 6;
Nnf = 4;
Naf = 4;

% Initialisation
ExtFai = zeros(Nf*size(Cn0,1),Nnf);            % Extrémités facettes en noeuds
ExtFaAri = zeros(Nf*size(Cn0,1),Naf);          % Extrémités facettes en arêtes
Cfi = [];                                     % Connectivités facettes

%Tab = [1 2 3 4; 5 8 7 6; 1 5 6 2; 2 6 7 3; 3 7 8 4; 4 8 5 1];
Tab = [1 2 3 4; 5 6 7 8; 2 3 7 6; 1 4 8 5; 1 5 6 2; 4 8 7 3];

for k = 1:1:Nf, ExtFai(k:Nf:Nf*size(Cn,1),:) = Cn(:,Tab(k,:)); end

%Tab = [1 2 3 4; 8 7 6 5; 9 5 10 1; 10 6 11 2; 11 7 12 3; 12 8 9 4];
Tab = [1 2 3 4; 5 6 7 8; 2 11 6 10; 4 12 8 9; 9 5 10 1; 12 7 11 3];

for k = 1:1:Nf, ExtFaAri(k:Nf:Nf*size(Cn,1),:) = abs(Ca(:,Tab(k,:))); end

for k = 1:1:Nf, Cfi = [Cfi, Nf*(1:size(Cn,1))'-(Nf-k)]; end

% Eliminer les doublons
[ExtFa,I,J] = unique(sort(ExtFai,2),'rows');
ExtFa = ExtFai(I,:);
Cf = J(Cfi);
Cf = int32(Cf);
if size(Cf,2) == 1, Cf = Cf'; end

if SaveOption == 1

[ExtFaAr,I,J] = unique(sort(ExtFaAri,2),'rows');
ExtFaAr = ExtFaAri(I,:);

% Signes des arêtes et coordonnées des facettes
TabCoorF = single(zeros(size(ExtFa,1),size(CoorN,2),size(ExtFa,2))); 
for k = 1:1:size(TabCoorF,3)
    TabCoorF(:,:,k) = single(CoorN(ExtFa(:,k),:));
end
CoorF = sum(TabCoorF,3)/size(TabCoorF,3);
clear TabCoorF
%
%CoorF = zeros(length(ExtFa),size(CoorN,2));
for k = 1:1:size(ExtFaAr,1)
    for ia = 1:1:size(ExtFaAr,2), if ExtAr(ExtFaAr(k,ia),1) ~= ExtFa(k,ia), ExtFaAr(k,ia) = -ExtFaAr(k,ia); end, end
    %CoorF(k,:) = sum(CoorN(ExtFa(k,:),:),1)/size(ExtFa,2); 
end
ExtFaAr = int32(ExtFaAr);
end

% -------------------------- Sous-domaines -------------------------%


NsdF = ones(size(ExtFa,1),1);
[xg,yg,zg] = deal(CoorV(:,1),CoorV(:,2),CoorV(:,3));

% Numéro des sous-domaines

Nsd = (length(l1)+1)*(length(l2)+1)*ones(size(Cn,1),1);    % domaine extrême

if nargin < 10
xi = -sum(l1)/2;
yi = -sum(l2)/2;
for k1 = 1:length(l1)
    for k2 = 1:length(l2)
        Nsd(xg>xi & xg<xi+l1(k1) & yg>yi & yg<yi+l2(k2)) = k1+k2*(length(l2)-1); 
        yi = yi+l2(k2);
    end
    xi = xi+l1(k1);
end
end



TypeElmt = 'Hexahedron';    % Type d'élément

Cn = int32(Cn);
Ca = int32(Ca);

if max(Nsd)<128, Nsd = int8(Nsd); else, Nsd = int16(Nsd); end
ExtAr = int32(ExtAr);
ExtFa = int32(ExtFa);


if SaveOption == 0
    Mesh = struct('Cn',Cn,'Ca',Ca,'Cf',Cf,'CoorN',CoorN,'CoorV',CoorV,...
              'Nsd',Nsd,'ExtFa',ExtFa,...
              'TypeElmt',TypeElmt,'NsdF',NsdF);
else
    Mesh = struct('Cn',Cn,'Ca',Ca,'Cf',Cf,'CoorN',CoorN,'CoorA',CoorA,'CoorF',CoorF,'CoorV',CoorV,...
              'Nsd',Nsd,'ExtAr',ExtAr,'ExtFa',ExtFa,'ExtFaAr',ExtFaAr,...
              'TypeElmt',TypeElmt,'NsdF',NsdF);
end
%
%
if ~isempty(R)
    Mesh.Nsd(:) = length(R)+1;
    r = sqrt(sum((Mesh.CoorV(:,1:2)).^2,2));
    for kr = 1:length(R)
        Mesh.Nsd(r<=R(kr)) = length(R)+1-kr;
    end
    %Mesh.R = R;
else
    if isempty(yv)
    if isempty(l1) && isempty(l2)
        Mesh.Nsd(:) = 1;
    elseif length(l1) == 1 && isempty(l2)
        Mesh.Nsd(:) = 2; Mesh.Nsd(abs(Mesh.CoorV(:,1))<l1/2) = 1;
    elseif length(l2) == 1 && isempty(l1)
        Mesh.Nsd(:) = 2; Mesh.Nsd(abs(Mesh.CoorV(:,2))<l2/2) = 1;
    elseif length(l1) == 1 && length(l2) == 1 
        Mesh.Nsd(:) = 2; Mesh.Nsd(abs(Mesh.CoorV(:,1))<l1/2 & abs(Mesh.CoorV(:,2))<l2/2) = 1;
    else
        rx = abs(x0(end-1:-1:floor(length(x0)/2)+1));
        ry = abs(y0(end-1:-1:floor(length(y0)/2)+1));
        %
        Mesh.Nsd(:) = length(find(rx>0))+1;
        for kr = 1:length(rx)
            if length(rx) == length(ry)
                P = abs(Mesh.CoorV(:,1))<=rx(kr) & abs(Mesh.CoorV(:,2))<=ry(kr);
            else
                P = abs(Mesh.CoorV(:,1))<=rx(kr);
            end
            Mesh.Nsd(P) = length(find(rx>0))+1-kr;
        end
    end
    end
        
end
%%
%%
if nargin <= 11
    if ~isempty(yv)
        if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(:) = 2; else, Mesh.Nsd(:) = max(NumSD)+1; end
        
        %Mesh.Nsd(:) = 2;
        % for k = 1:size(xv,1)
        % in = inpolygon(Mesh.CoorV(:,1),Mesh.CoorV(:,2),xv(k,:),yv(k,:));
        % if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in) = 1; else, Mesh.Nsd(in) = NumSD(k); end
        % %Mesh.Nsd(in) = 1;
        % end
        % in = cell(size(xv,1),1);
        % [xf,yf] = deal(Mesh.CoorV(:,1),Mesh.CoorV(:,2));
        % parfor k = 1:size(xv,1)
        %     in{k} = inpolygon(xf,yf,xv(k,:),yv(k,:));
        % end
        %tic
        if ~isempty(yv) && npz==2,
            [xv,yv]=deal(single(xv),single(yv));
            Xg = sum(xv,2)/size(xv,2);
            Yg = sum(yv,2)/size(yv,2);
            
            [xf,yf] = deal(single(Mesh.CoorV(:,1)),single(Mesh.CoorV(:,2)));
            
            P = int32(1:length(xf));
            for k = 1:size(xv,1)
                P0 = (abs(xf-Xg(k))<1.01*dx(k) & abs(yf-Yg(k))<1.01*dy(k));
                Pn = P(P0);
                ik = inpolygon(xf(Pn),yf(Pn),xv(k,:),yv(k,:));
                in = Pn(ik);
                if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in) = 1; else, Mesh.Nsd(in) = NumSD(k); end
            end
        else
            for k = 1:size(xv,1)
                in = inpolygon(Mesh.CoorV(:,1),Mesh.CoorV(:,2),xv(k,:),yv(k,:));
                if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in) = 1; else, Mesh.Nsd(in) = NumSD(k); end
            end

        end
%toc
        
        % for k = 1:size(xv,1)
        %     if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in{k}) = 1; else, Mesh.Nsd(in{k}) = NumSD(k); end
        % end

    end
    Mesh.xv = xv;
    Mesh.yv = yv;

end

%%
% Recherche des facettes extérieures
Pf = zeros(length(Mesh.NsdF),1);%zeros(length(Mesh.CoorF),1);
NsdF = zeros(length(Mesh.NsdF),1);%zeros(length(Mesh.CoorF),1);
%tic
for ie = 1:size(Mesh.Cf,1)
    for k = 1:size(Mesh.Cf,2)
        if Pf(Mesh.Cf(ie,k)) == 0 
            Pf(Mesh.Cf(ie,k)) = 1;
            NsdF(Mesh.Cf(ie,k)) = Mesh.Nsd(ie);
        else
            Pf(Mesh.Cf(ie,k)) = Pf(Mesh.Cf(ie,k))+1;
            NsdF(Mesh.Cf(ie,k)) = 0;
        end
    end
end
%toc

%
Mesh.NsdF = int16(NsdF);
Mesh.TabP = int16(Pf);
%
% Volume des éléments
Mesh.Vol = zeros(size(Mesh.CoorV,1),1);

d21 = sqrt(sum((Mesh.CoorN(Mesh.Cn(:,2),:)-Mesh.CoorN(Mesh.Cn(:,1),:)).^2,2));
d41 = sqrt(sum((Mesh.CoorN(Mesh.Cn(:,4),:)-Mesh.CoorN(Mesh.Cn(:,1),:)).^2,2));
d51 = sqrt(sum((Mesh.CoorN(Mesh.Cn(:,5),:)-Mesh.CoorN(Mesh.Cn(:,1),:)).^2,2));

Mesh.Vol = d21.*d41.*d51;

Mesh = NumFacet(Mesh);

end 

%% ------------------------------------------------------------------------%


function Mesh1 = NumFacet(Mesh)

for k = 1:max(Mesh.Nsd), s(k).TabP = int16(zeros(length(Mesh.ExtFa),1)); end


for kn = 1:max(Mesh.Nsd),

    Pe = find(Mesh.Nsd == kn);
    %Pfe = unique(Mesh.Cf(Pe,:));
    %TabPf = zeros(length(Pfe),1);
    
    % for ke = 1:length(Pe),
    %     for k = 1:size(Mesh.Cf,2),
    %         s(kn).TabP(Mesh.Cf(Pe(ke),k)) = s(kn).TabP(Mesh.Cf(Pe(ke),k))+1;
    %     end
    % end
    % for ke = 1:length(Pe)
    %       s(kn).TabP(Mesh.Cf(Pe(ke),:)) = s(kn).TabP(Mesh.Cf(Pe(ke),:))+1;
    % end
    s(kn).TabP(Mesh.Cf(Pe(:),:)) = s(kn).TabP(Mesh.Cf(Pe(:),:))+1;

    s(kn).TabP(s(kn).TabP ~= 1) = int16(0);
    s(kn).TabP(s(kn).TabP == 1) = int16(kn);
    
    Mesh.TabNsdF{kn} = int16(s(kn).TabP);

%    disp(kn)

end

Mesh1 = Mesh;

end

%% ------------------------------------------------------------------------%
function [Ca,ExtAr,CoorA,CoorF] = RechercheArete(Cn,CoorN)

% RechercheArete
%       Recherche des informations sur les arêtes
%
% Syntaxe
%   [Ca,ExtAr,CoorA,CoorF] = RechercheArete(Cn,CoorN);
%
% Description
%   Cn : connectivité aux noeuds 
%   CoorN : coordonnées aux noeuds
% 
%   Ca : connectivités arêtes
%   ExtAr : tableau des noeuds extrimités des arêtes 
%   CoorA : coordonnées des arêtes
%   CoorF : coordonnées des facettes

%   Date de la dernière version : 18 janvier 2019
%   Auteur : Mondher Besbes (LCF / CNRS / IOGS)


% Recherche des facettes
% ----------------------

% Initialisation
ExtAri = zeros(numel(Cn),2);                % Extrémités arêtes
CoorF = zeros(length(Cn),size(CoorN,2));    % Coordonnées des facettes
Cai = [];                                   % Connectivités arêtes

for k = 1:1:size(Cn,2), ExtAri(k:size(Cn,2):numel(Cn),:) = Cn(:,[k k*(k<size(Cn,2))+1]); end
for k = 1:1:size(Cn,2), Cai = [Cai, size(Cn,2)*(1:size(Cn,1))'-(size(Cn,2)-k)]; end

% Eliminer les doublons
[ExtAr,I,J] = unique(sort(ExtAri,2),'rows');
ExtAr = ExtAri(I,:);
Ca(1:size(Cn,1),:) = J(Cai);

% Coordonnées des arêtes
%CoorA = zeros(length(ExtAr),size(CoorN,2));
%for k = 1:1:length(ExtAr), CoorA(k,:) = sum(CoorN(ExtAr(k,:),:),1)/size(ExtAr,2); end
TabCoorA = zeros(length(ExtAr),size(CoorN,2),size(ExtAr,2)); 
for k = 1:1:size(TabCoorA,3)
    TabCoorA(:,:,k) = CoorN(ExtAr(:,k),:);
end
CoorA = sum(TabCoorA,3)/size(TabCoorA,3);
                   
% Signes des arêtes et coordonnées des facettes
%CoorF = zeros(size(Cn,1),size(CoorN,2));
TabCoorF = zeros(size(Cn,1),size(CoorN,2),size(Cn,2)); 
for k = 1:1:size(TabCoorF,3)
    TabCoorF(:,:,k) = CoorN(Cn(:,k),:);
end
CoorF = sum(TabCoorF,3)/size(TabCoorF,3);

for ie = 1:size(Cn,1),
    for ia = 1:1:size(Ca,2), if ExtAr(abs(Ca(ie,ia)),1) ~= Cn(ie,ia), Ca(ie,ia) = -Ca(ie,ia); end, end
    %CoorF(ie,:) = sum(CoorN(Cn(ie,:),:),1)/size(Cn,2);
end

end

%% ----------------------------------------------------------------------
function [x0,y0,dx,dy] = OptimMesh(x0,y0,xv,yv,NumSD)


% Mesh 2D
[xv,yv]=deal(single(xv),single(yv));
[x0,y0]=deal(single(x0),single(y0));

[x,y] = meshgrid(x0,y0);

Cn = int32(zeros((length(x0)-1)*(length(y0)-1),4));%
%Cn = [];
for i = 1:length(y0)-1
    i0 = (1:length(x0)-1) + (i-1)*length(x0);
    i1 = i0 + length(x0);
    %Cn = [Cn; [i0' i0'+1 i1'+1 i1']];
    Cn((1:length(x0)-1)+(i-1)*(length(x0)-1),:) = int32([i0' i0'+1 i1'+1 i1']);
end

    
%CoorN = [];
CoorN = single(nan((length(x0))*(length(y0)),2));
j = 1:length(x0);
for i = 1:length(y0)
    %CoorN = [CoorN; [x(i,j)' y(i,j)']];
    CoorN(j+(i-1)*length(x0),:) = [x(i,j)' y(i,j)'];
end

% Numéro des sous-domaines

Nsd = int16(ones(size(Cn,1),1));    % domaine extrême

% Recherche des numéros des arêtes ...
% [Ca,ExtAr,CoorA,CoorF] = RechercheArete(Cn,CoorN);
% TypeElmt = 'Quadrangle';
% 
% % Structure Mesh
% Mesh = struct('Cn',Cn,'Ca',Ca,'CoorN',CoorN,'CoorA',CoorA,'CoorF',CoorF,...
%               'Nsd',Nsd,'ExtAr',ExtAr,'TypeElmt',TypeElmt);

TabCoorF = single(zeros(size(Cn,1),size(CoorN,2),size(Cn,2))); 
for k = 1:1:size(TabCoorF,3)
    TabCoorF(:,:,k) = CoorN(Cn(:,k),:);
end
CoorF = sum(TabCoorF,3)/size(TabCoorF,3);

clear Cn CoorN

%if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(:) = 2; else, Mesh.Nsd(:) = max(NumSD)+1; end
if isempty(NumSD) || length(NumSD) == 1, Nsd(:) = 2; else, Nsd(:) = max(NumSD)+1; end

%
% for k = 1:size(xv,1)
%     in = inpolygon(Mesh.CoorF(:,1),Mesh.CoorF(:,2),xv(k,:),yv(k,:));
%     if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in) = 1; else, Mesh.Nsd(in) = NumSD(k); end
% end

%tic
Xg = sum(xv,2)/size(xv,2);
Yg = sum(yv,2)/size(yv,2);

Xv = (xv-repmat(Xg,1,size(xv,2)));
Yv = (yv-repmat(Yg,1,size(xv,2)));
dx = max(abs(Xv),[],2);
dy = max(abs(Yv),[],2);

%in = cell(size(xv,1),1);
[xf,yf] = deal(CoorF(:,1),CoorF(:,2));

P = int32(1:length(xf))';
for k = 1:size(xv,1)
    P0 = (abs(xf-Xg(k))<1.01*dx(k) & abs(yf-Yg(k))<1.01*dy(k));
    Pn = P(P0);
    ik = inpolygon(xf(Pn),yf(Pn),xv(k,:),yv(k,:));
    in = Pn(ik);
    if isempty(NumSD) || length(NumSD) == 1, Nsd(in) = 1; else, Nsd(in) = NumSD(k); end
end

%toc
clear CoorF xf yf Xv Yv


% for k = 1:size(xv,1)
%     %if isempty(NumSD) || length(NumSD) == 1, Mesh.Nsd(in{k}) = 1; else, Mesh.Nsd(in{k}) = NumSD(k); end
%     if isempty(NumSD) || length(NumSD) == 1, Nsd(in{k}) = 1; else, Nsd(in{k}) = NumSD(k); end
% end


%figure, VisuMesh(Mesh)
%u = reshape(Mesh.Nsd,length(x0)-1,length(y0)-1);
u = reshape(Nsd,length(x0)-1,length(y0)-1);

% Reduce Mesh
[mx,my] = size(u);
%
x = x0(mx+1);
%uu = u(mx,:);
uu = nan(size(u));
uu(mx,:) = u(mx,:);
for ii = mx-1:-1:1
    if ~all(u(ii,:) == u(ii+1,:))
        %uu = [u(ii,:);uu];
        uu(ii,:) = u(ii,:);
        x = [x0(ii+1);x];
    end
end

P = ~isnan(uu(:,1));
uu = uu(P,:);

u = uu(:,my);
y = y0(my+1);
for jj = my-1:-1:1
    if ~all(uu(:,jj)==uu(:,jj+1))
%        u = [uu(:,jj),u];
        y = [y0(jj+1);y];
    end
end

x0 = [x0(1);x]';
y0 = [y0(1);y]';

end




