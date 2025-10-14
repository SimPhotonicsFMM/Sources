function [E,H] = CalculFieldFMM(Data,Mesh,Phys,Sb,MatS,Sh,X,Y,Z)

% CalculFieldMMF
%   Electromagnetic field computation at point, line or grid of points
%
% Syntax
%   [E,H] = CalculFieldFMM(Model,X,Y,Z,pol); % structure 3D
%   [E,H] = CalculFieldFMM(Model,X,Y);       % structure 2D
%   [E,H] = CalculFieldFMM(Data,Mesh,Phys,Sb,MatS,Sh,X,Y,Z); % structure 3D
%   [E,H] = CalculFieldFMM(Data,Mesh,Phys,Sb,MatS,Sh,X,Y);   % structure 2D
%
% Description
%   Model : Model structure from "spectrum" function
%   Data : Data of the studied problem (wavelength, angle of incidence...)
%   Mesh : Mesh Structure 
%   Phys : Physical properties (see CaractMat.m)
%   Sb   : S-matrix of Lower middle 
%   Sh   : S-matrix of upper middle
%   MatS : S-matrix of different layers
%   pol  : Polarization 0: TE , 2:TM (by default the value in Data)
%   X,Y,Z : Array or Grid of computation points
%
%   E : Electric field
%   H : Magnetic field
%
% Example
%   lambda = 0.85; theta = 56*pi/180; inc = -1;
%   dx = 0.2; dy = .2; geom = [.015 0.04];
%   index = {1.33 , {IndexVal('Au') 1.33} ,IndexVal('Au'), 1.7};
%   s = Spectrum(index,geom,lambda,theta,inc,...
%           'dx',dx,'dy',dy,'lix',{dx/2 0},'liy',{dy/2 0},'mx',10,'my',10);
%   [x,y,z] = ndgrid(linspace(-dx/2,dx/2,100),linspace(-dy/2,dy/2,101),0);
%   [VectE,VectH] = CalculFieldFMM(s,x,y,z); 

% Date of the latest version : 02 May 2024
% Author : Mondher Besbes (LCF / CNRS / IOGS)


if nargin == 4 || nargin == 5, Z = Sb; end
if nargin <= 5
    s = Data;
    if nargin == 5, pol = MatS; end
    if nargin == 2
        [x,y] = deal(Mesh(:,1),Mesh(:,2));
        z = Mesh(:,end);
    else
        [X,Y] = deal(Mesh,Phys);
    end
    [Data,Mesh,Phys,Sb,MatS,Sh] = deal(s.Data,s.Mesh,s.Phys,s.Sb,s.MatS,s.Sh);
    if ~isfield(s,'Nper'), s.Nper = 1; end
    Mesh = repmat(Mesh,1,s.Nper);
    for k = 2:s.Nper
        if size(Mesh(1).CoorN,2) == 3
            Mesh(k) = MoveMesh(utilMesh(Mesh(k)),[0 0 max(Mesh(k-1).CoorN(:,end))]);
        else
            Mesh(k) = MoveMesh(utilMesh(Mesh(k)),[0 max(Mesh(k-1).CoorN(:,end))]);
        end
    end
    MatS = repmat(MatS,s.Nper,1);
    Phys = repmat(Phys,1,s.Nper);
    Data = repmat(Data,1,s.Nper);
end
if nargin == 5
    for k = 1:length(Data)
     Data(k).TypePol = pol;
     Phys(k).TypePol = pol;
    end
end
%
if nargin == 9 || nargin == 4 || nargin == 5
    if isvector(X) && isvector(Y) && isvector(Z)
        [x,y,z] = ndgrid(X,Y,Z);
    else
        [x,y,z] = deal(X,Y,Z);
    end
elseif nargin == 8 || nargin == 3
    if isvector(X) && isvector(Y)
        [x,y] = ndgrid(X,Y);
    else
        [x,y] = deal(X,Y);
    end
    z = y;
end
%
if sum(Data(1).Sym) ~= 4
    disp('For Field Calculation, symmetry is not considered');
end
%
[E,H] = deal(nan(length(x(:)),6));

a0 = max(Mesh(1).CoorN(:,1))-min(Mesh(1).CoorN(:,1));
b0 = max(Mesh(1).CoorN(:,2))-min(Mesh(1).CoorN(:,2));

% Milieu bas
if min(z(:))<min(Mesh(1).CoorN(:,end))
    %
    h0 = min(Mesh(1).CoorN(:,end))-min(z(:));
    if abs(Data(1).Lambda0/4/Data(1).nb-h0)<10*eps, h0 = h0+h0/1e6; end
    %
    if size(Mesh(1).CoorN,2) == 3
        Mesh0 = MeshLayer(a0,[],b0,[],h0,2,2,2);
        %Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 0 min(z(:))+min(Mesh(1).CoorN(:,end))]);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 0 min(Mesh(1).CoorN(:,end))-h0]);
    elseif size(Mesh(1).CoorN,2) == 2
        Mesh0 = MeshLayer(a0,[],h0,2,2);
        %Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 min(z(:))+min(Mesh(1).CoorN(:,end))]);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 min(Mesh(1).CoorN(:,end))-h0]);
    end
    %
    Data0 = SetData(Data(1),'Indice',Data(1).nb);
    Phys0 = CaractMat(Mesh0,Data0,Phys(1).Kx,Phys(1).Ky);
    S0 = CalculMatS(Data0,Mesh0,Phys0); % milieu bas 
    if isfield(Mesh(1),'xv'), Mesh0.xv = []; Mesh0.yv = []; end
    if isfield(Mesh(1),'zv'), Mesh0.zv = []; end
    Mesh = [Mesh0 Mesh];
    MatS = [S0; MatS];
    Phys = [Phys0 Phys];
    Data = [Data0 Data];
end

% Milieu haut
if max(z(:))>max(Mesh(end).CoorN(:,end))
    %
    h0 = max(z(:))-max(Mesh(end).CoorN(:,end));
    %
    if size(Mesh(1).CoorN,2) == 3
        Mesh0 = MeshLayer(a0,[],b0,[],h0,2,2,2);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 0 max(Mesh(end).CoorN(:,end))]);
    elseif size(Mesh(1).CoorN,2) == 2
        Mesh0 = MeshLayer(a0,[],h0,2,2);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 max(Mesh(end).CoorN(:,end))]);
    end
    %
    Data0 = SetData(Data(1),'Indice',Data(1).nh);
    Phys0 = CaractMat(Mesh0,Data0,Phys(1).Kx,Phys(1).Ky);
    S0 = CalculMatS(Data0,Mesh0,Phys0); % milieu haut 
    if isfield(Mesh(1),'xv'), Mesh0.xv = []; Mesh0.yv = []; end
    if isfield(Mesh(1),'zv'), Mesh0.zv = []; end
    Mesh = [Mesh Mesh0];
    MatS = [MatS; S0];
    Phys = [Phys Phys0];
    Data = [Data Data0];
end

%tic
for k = 1:length(Mesh)
    % Grille de calcul
    Pn = find(z(:)>=min(Mesh(k).CoorN(:,end)) & z(:)<=max(Mesh(k).CoorN(:,end))+eps/2);
    %
    if ~isempty(Pn)
        x0 = x(Pn); y0 = y(Pn); z0 = z(Pn);
        %
        if abs(min(z0)-min(Mesh(k).CoorN(:,end)))>10*eps
            x0 = [x0(:);min(Mesh(k).CoorN(:,1))];
            y0 = [y0(:);min(Mesh(k).CoorN(:,2))];
            z0 = [z0(:);min(Mesh(k).CoorN(:,end))];
        end
        if abs(max(z0(:))-max(Mesh(k).CoorN(:,end)))>10*eps
            x0 = [x0(:);max(Mesh(k).CoorN(:,1))];
            y0 = [y0(:);max(Mesh(k).CoorN(:,2))];
            z0 = [z0(:);max(Mesh(k).CoorN(:,end))];
        end
        %
        % Mise à jour de Sb1 et Sh1
        Sb1 = Sb;
        for k1 = 1:k-1, Sb1 = ProdMatS(Sb1,MatS(k1,:)); end
        %
        Sh1 = Sh;
        for k1 = length(Mesh):-1:k+1, Sh1 = ProdMatS(MatS(k1,:),Sh1); end

        % Calcul du champ par cellule
%         [VectE,VectH] = FieldFMM(Data(k),Mesh(k),Phys(k),Sb1,MatS(k,:),Sh1,[x0(:) y0(:) z0(:)]);
%         % 
%         E(Pn,1:end) = VectE{1}(1:length(Pn),:);
%         H(Pn,1:end) = VectH{1}(1:length(Pn),:);
        %
        if size(Mesh(1).CoorN,2) == 2
            [VectE,VectH] = FieldMMF2D(Data(k),Mesh(k),Phys(k),Sb1,MatS(k,:),Sh1,[x0(:) y0(:) z0(:)]);
        elseif size(Mesh(1).CoorN,2) == 3
            [VectE,VectH] = FieldMMF3D(Data(k),Mesh(k),Phys(k),Sb1,MatS(k,:),Sh1,[x0(:) y0(:) z0(:)]);
        end
        %
        E(Pn,1:end) = VectE(1:length(Pn),:);
        H(Pn,1:end) = VectH(1:length(Pn),:);

    end
end
%toc

end

% ------------------------------------------------------------------------%

function [VectE,VectH,Vect] = FieldFMM(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor)

% FieldFMM
%   Calcul du champ électromagnétique aux noeuds du maillage (une grille)
%
% Syntaxe
%   [VectE,VectH,Vect] = FieldFMM(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor);
%   [VectE,VectH,Vect] = FieldFMM(Data,Mesh,Phys,Sb,MatS,Sh);
%   [VectE,VectH,Vect] = FieldFMM(Data,Mesh,Phys);
%
% Description
%   Mesh : structure données du maillage (CoorN, Cn, CoorA, Ca, ...,Nsd )
%   Phys : structure données matériaux (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%   Sb : Matrice S du milieu bas
%   Sh : Matrice S du milieu haut
%   MatS : Matrice S des différentes tranches
%
%   VectE : Vecteur champ électrique
%   VectH : Vecteur champ magnétique
%   Vect :   TE : Hx, Hy, Ez
%                TM : Ex, Ey, Hz
%
% Example
%   Same example in the function "CalculMatS"
%   Data = SetData(Data,'TypePol',2);   % 0:TE  2:TM
%   Phys = CaractMat(Mesh,Data);
%   [x,y,z] = ndgrid(linspace(-dx/2,dx/2,100),linspace(-dy/2,dy/2,101),h/2); % (xoy) 
%   [VectE,VectH] = FieldMMF(Data,Mesh,Phys,Sb,MatS,Sh,[x(:) y(:) z(:)]); 

% Date de la dernière version : 24 août 2021
% Auteur : Mondher Besbes (LCF / CNRS / IOGS)

if nargin < 5
    Sb = CalculMatS(Data,Mesh,Phys,-1); % milieu bas
    Sh = CalculMatS(Data,Mesh,Phys,+1); % milieu haut
end
%
if nargin < 4
    MatS = CalculMatS(Data,Mesh,Phys); % Matrices S des différentes couches
end

if nargin < 7, TabCoor = []; end

Sb1 = Sb;
 [VectE,VectH,Vect] = deal(cell(size(MatS,1),1));
for kc = 1:size(MatS,1)-1
    S1 = ProdMatS(MatS(kc+1:end,:));
    Sh1 = ProdMatS(S1,Sh);
    if size(Mesh(1).CoorN,2) == 2
        [VectE{kc},VectH{kc},Vect{kc}] = FieldMMF2D(Data(kc),Mesh(kc),Phys(kc),Sb1,MatS(kc,:),Sh1,TabCoor);
    elseif size(Mesh(1).CoorN,2) == 3
        [VectE{kc},VectH{kc},Vect{kc}] = FieldMMF3D(Data(kc),Mesh(kc),Phys(kc),Sb1,MatS(kc,:),Sh1,TabCoor);
    end
    S1 = ProdMatS(MatS(1:kc,:));
    Sb1 = ProdMatS(Sb,S1);
end
if size(Mesh(1).CoorN,2) == 2
    [VectE{end},VectH{end},Vect{end}] = FieldMMF2D(Data(end),Mesh(end),Phys(end),Sb1,MatS(end,:),Sh,TabCoor);
elseif size(Mesh(1).CoorN,2) == 3
    [VectE{end},VectH{end},Vect{end}] = FieldMMF3D(Data(end),Mesh(end),Phys(end),Sb1,MatS(end,:),Sh,TabCoor);
end
%
end

%% ------------------------------------------------------------------------
% Structure 2D
% -------------------------------------------------------------------------

function [VectE,VectH,Vect] = FieldMMF2D(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor)

% FieldMMF2D
%   Calcul du champ électromagnétique aux noeuds du maillage (une grille)
%
% Syntaxe
%   [VectE,VectH,Vect] = FieldMMF2D(Data,Mesh,Phys,Sb,MatS,Sh);
%
% Description
%   Mesh : structure données du maillage (CoorN, Cn, CoorA, Ca, ...,Nsd )
%   Phys : structure données matériaux (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%
%   VectE : vecteur champ électrique
%   VectH : vecteur champ magnétique
%              TE : Hx, Hy, Ez
%              TM : Ex, Ey, Hz

% Date de la dernière version : 12 septembre 2017
%                               21 octobre 2021 (Apod =1 si mx = 0)
%                               04 mars 2022 (Calcul sur une grille)
% Auteur : Mondher Besbes (LCF / CNRS / IOGS)

%
if isempty(TabCoor)
    [x,y] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2));
    y0 = y;
else
    Pn = TabCoor(:,2)>=min(Mesh.CoorN(:,2)) & TabCoor(:,2)<=max(Mesh.CoorN(:,2));
    [x,y] = deal(TabCoor(Pn,1),TabCoor(Pn,2));
    y0 = y;
end

hc = max(Mesh.CoorN(:,2))-min(Mesh.CoorN(:,2));%Data.hc; % hauteur de la couche 
ymin = min(Mesh.CoorN(:,2));
if abs(hc-(max(y0)-min(y0)))>10*eps 
    error('Hauteur de la grille est différente de la hauteur de la couche')
end
%
% hc = max(Mesh.CoorN(:,2))-min(Mesh.CoorN(:,2));%Data.hc; % hauteur de la couche 
% ymin = min(Mesh.CoorN(:,2));
m = size(Sb{3},1);
mx = (m-1)/2;

Dx = max(Mesh.CoorN(:,1))-min(Mesh.CoorN(:,1));

Beta0 = Phys.Kx;
BetaX = Beta0 + (-mx:1:mx)*2*pi/Dx;

% Calcul du champ

[Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});

[Ib,Ih] = deal(zeros(m,1));
if Data.ChampInc == -1, Ib(mx+1) = 1; elseif Data.ChampInc == +1, Ih(mx+1) = 1; end

MatSh1 = ProdMatS(MatS,Sh);
Sh1 = MatSh1{1};
[nh1,nh2,nh3,nh4] = deal(m,length(Pui),length(Pu),m);
[Sh11,Sh12,Sh21,Sh22] = deal(Sh1(1:nh3,1:nh1),Sh1(1:nh3,nh1+1:end),...
                             Sh1(nh3+1:end,1:nh1),Sh1(nh3+1:end,nh1+1:end));
Sb1 = Sb{1};
[nb1,nb2,nb3,nb4] = deal(length(Pdi),m,m,length(Pd));
[Sb11,Sb12,Sb21,Sb22] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end),...
                             Sb1(nb3+1:end,1:nb1),Sb1(nb3+1:end,nb1+1:end));
                         
SIh = Sh22*Ih(Pui-m);
if size(SIh,2) == 0, SIh = zeros(size(SIh,1),1); end
SIb = Sb11*Ib(Pdi);
if size(SIb,2) == 0, SIb = zeros(size(SIb,1),1); end

EHb = [-Sh21 eye(m) ; eye(m) -Sb12]\[SIh ; SIb];
%

[Ex,Ey,Hz] = deal(zeros(1,length(x)));

%

[Em0,Hm0] = deal(EHb(1:m),EHb(1+m:end));
P = MatS{3};
Q = MatS{4};
InvP = inv(P);
InvQ = inv(Q);
%
Vp = MatS{6};
TpzMatQy = MatS{7};
InvTpzMatQy = inv(TpzMatQy);

MatQP = [Q , Q; P -P];
InvMatQP = 0.5*[InvQ , InvP; InvQ -InvP];


% try
%     dP = decomposition(P,'lu');
%     dQ = decomposition(Q,'lu');
% catch
%     dP = P;
%     dQ = Q;
% end


%
[Em,Hm] = deal(Em0,Hm0);
[MatSb1,MatSh1] = deal(MatS);
% Apodisation : fenêtre de Hamming
if mx == 0
    TfApod = 1;
else
    TfApod = (0.54-0.46*cos((BetaX(:)+min(BetaX))*2*pi/(max(BetaX)-min(BetaX))));
end
%
C = InvMatQP*EHb;%(EHb.*repmat(TfApod,2,1));
%
CondM = rcond(InvMatQP);
MatQyBetaX = InvTpzMatQy*diag(BetaX);
%
while ~isempty(y)
    Pn0 = abs(y0 - min(y)) <= max(abs(Mesh.CoorN(:)))/1e8;
    % Apodisation : fenêtre de Hamming
    Em = Em.*TfApod;
    Hm = Hm.*TfApod;
    %
    Phase = exp(1i*x(Pn0)*BetaX);
    if isfield(Data,'DofType') && strcmp(Data.DofType,'Edge')
        Eym = -1i/Phys.K0*MatQyBetaX*Hm;
        Ex(Pn0) = Phase*Em;
        Hz(Pn0) = Phase*Hm;
        Ey(Pn0) = Phase*Eym;
    else
        Hym = -1i/Phys.K0*MatQyBetaX*Em;
        Ez(Pn0) = Phase*Em;
        Hx(Pn0) = Phase*Hm;
        Hy(Pn0) = Phase*Hym;
    end
    %
    Pn = y > min(y);
    if isempty(y(Pn)), break; end
    hy = min(y(Pn))-ymin;
    %
    ExpVp = exp([Vp ;-Vp]*hy);
    if 1/norm(ExpVp(end/2+1:end))*CondM < eps
    %if rcond(repmat([ones(size(Vp));exp(-Vp*hy)],1,2*length(Vp)).*InvMatQP) < eps/1e3 %max(abs(exp(-Vp*hy))) > 1e3/eps, % % hy>=Data.Lambda0/4/n0 %  %hy>=Data.Lambda0/4 % 
        hy0 = hy + ymin - min(y0);
        InvCosH = repmat(reshape(1./cosh(Vp*hy0),1,length(Vp)),length(Vp),1);
        TanH = repmat(reshape(tanh(Vp*hy0),1,length(Vp)),length(Vp),1);
        MatSb1{1} = [(Q.*InvCosH)*InvQ , (Q.*TanH)*InvP; -(P.*TanH)*InvQ , (P.*InvCosH)*InvP];
        %
        InvCosH = repmat(reshape(1./cosh(Vp*(hc-hy0)),1,length(Vp)),length(Vp),1);
        TanH = repmat(reshape(tanh(Vp*(hc-hy0)),1,length(Vp)),length(Vp),1);
        MatSh1{1} = [(Q.*InvCosH)*InvQ , (Q.*TanH)*InvP; -(P.*TanH)*InvQ , (P.*InvCosH)*InvP];
        %
        [MatSb1{2},MatSh1{2}] = deal(MatS{2});
        MatSh1 = ProdMatS(MatSh1,Sh);
        Sh1 = MatSh1{1};
        [Sh21,Sh22] = deal(Sh1(nh3+1:end,1:nh1),Sh1(nh3+1:end,nh1+1:end));
        MatSb1 = ProdMatS(Sb,MatSb1);
        Sb1 = MatSb1{1};
        [Sb11,Sb12] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end));
    
        SIh = Sh22*Ih(Pui-m);
        if size(SIh,2) == 0, SIh = zeros(size(SIh,1),1); end
        SIb = Sb11*Ib(Pdi);
        if size(SIb,2) == 0, SIb = zeros(size(SIb,1),1); end
    
        EHm = [-Sh21 eye(m) ; eye(m) -Sb12]\[SIh ; SIb];
        %
        C = InvMatQP*EHm;
        ymin = min(y(Pn));

    else
    %
        %ExpVp = exp([Vp ;-Vp]*hy);
        EHm = MatQP*(ExpVp.*C);
    end
    [Em,Hm] = deal(EHm(1:m),EHm(1+m:end));
    %
    y = y(Pn);
end


[VectE,VectH] = deal(zeros(length(x),3));

% attention au clonage
if isfield(Data,'DofType') && strcmp(Data.DofType,'Edge')
    if Phys.TypePol == 0,       % E
        VectE(:,3) = Hz(:)/(1i*sqrt(Phys.Z0));
        VectH(:,1:2) = [Ex(:) Ey(:)]*sqrt(Phys.Z0);
        Vect = [VectH(:,1:2) VectE(:,3)];
    elseif Phys.TypePol == 2,   % H
        VectE(:,1:2) = [Ex(:) Ey(:)]*sqrt(Phys.Z0);
        VectH(:,3) = Hz(:)/(1i*sqrt(Phys.Z0));
        Vect = [VectE(:,1:2) VectH(:,3)];
    end
else
    if Phys.TypePol == 0,       % E
        VectE(:,3) = Ez(:)*sqrt(Phys.Z0);
        VectH(:,1:2) = [Hx(:) Hy(:)]/(1i*sqrt(Phys.Z0));
        Vect = [VectH(:,1:2) VectE(:,3)];
    elseif Phys.TypePol == 2,   % H
        VectE(:,1:2) = [Hx(:) Hy(:)]/(1i*sqrt(Phys.Z0));
        VectH(:,3) = Ez(:)*sqrt(Phys.Z0);
        Vect = [VectE(:,1:2) VectH(:,3)];
    end
end


end


%% -------------------------------------------------------------------------
% Structure 3D
% -------------------------------------------------------------------------

function [VectE,VectH,Vect] = FieldMMF3D(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor)

if isempty(TabCoor)
    [x,y,z] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2),Mesh.CoorN(:,3));
    z0 = z;
else
    Pn = TabCoor(:,3)>=min(Mesh.CoorN(:,3)) & TabCoor(:,3)<=max(Mesh.CoorN(:,3))+eps/2;
    [x,y,z] = deal(TabCoor(Pn,1),TabCoor(Pn,2),TabCoor(Pn,3));
    z0 = z;
end

hc = max(Mesh.CoorN(:,3))-min(Mesh.CoorN(:,3));%Data.hc; % hauteur de la couche 
zmin = min(Mesh.CoorN(:,3));
if abs(hc-(max(z0)-min(z0)))>10*eps 
    error('Hauteur de la grille est différente de la hauteur de la couche')
end
%hc = ;%Data.hc; % hauteur de la couche 
%zmin = min(z0);      % Attention si les limites ne coïncident pas les min max de Mesh.CoorN(:,3)


[mx,my] = ndgrid(1:2*Data.mx+1,1:2*Data.my+1);
mx = mx(:)';
my = my(:)';
%
mx0 = Data.mx+1;
my0 = Data.my+1;
%
Dx = max(Mesh.CoorN(:,1))-min(Mesh.CoorN(:,1)); % Période en x
Dy = max(Mesh.CoorN(:,2))-min(Mesh.CoorN(:,2)); % Période en y
  

BetaX0 = Phys.Kx;
BetaY0 = Phys.Ky;
%
BetaX = BetaX0 + 2*pi/Dx*(mx-mx0);
BetaY = BetaY0 + 2*pi/Dy*(my-my0);

% Intensité incidente
m = size(Sb{3},1); 
[Ib,Ih] = deal(zeros(2*m,2)); % % dim 2 A vérifier

mx = m/2;
[Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});

if Data.ChampInc == +1
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY) 
            Ih(Sh{7}(1),1) = 1;
            Ih(Sh{7}(2),2) = 1;
        else
            if sum(Data.Sym) ~= 4
                [Ib,Ih] = deal(zeros(m,1));
                Ih(Pui) = 1;
            else
                Ih((mx-1)/2+1+3*mx,1) = 1; 
                Ih((mx-1)/2+1+2*mx,2) = 1;
            end
        end
else
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY), 
            Ib(Sb{7}(1),1) = 1; 
            Ib(Sb{7}(2),1) = 1;
        else
            Ib((mx-1)/2+1+mx,1) = 1; %Pdi = [(length(mx)-1)/2+1+length(mx) (length(mx)-1)/2+1];
            Ib((mx-1)/2+1,2) = 1;
        end
end
%

%MatSh1 = ProdMatS(MatS,Sh);
MatSh1 = Sh; for k1 = size(MatS,1):-1:1, MatSh1 = ProdMatS(MatS(k1,:),MatSh1); end

Sh1 = MatSh1{1};
[nh1,nh2,nh3,nh4] = deal(m,length(Pui),length(Pu),m);
[Sh11,Sh12,Sh21,Sh22] = deal(Sh1(1:nh3,1:nh1),Sh1(1:nh3,nh1+1:end),...
                             Sh1(nh3+1:end,1:nh1),Sh1(nh3+1:end,nh1+1:end));
Sb1 = Sb{1};
[nb1,nb2,nb3,nb4] = deal(length(Pdi),m,m,length(Pd));
[Sb11,Sb12,Sb21,Sb22] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end),...
                             Sb1(nb3+1:end,1:nb1),Sb1(nb3+1:end,nb1+1:end));
                         
SIh = Sh22*Ih(Pui,:);  % A vérifier
if size(SIh,2) == 0, SIh = zeros(size(SIh,1),size(Ih,2)); end
SIb = Sb11*Ib(Pdi,:);
if size(SIb,2) == 0, SIb = zeros(size(SIb,1),size(Ih,2)); end

%EHb = full([-Sh21 eye(m) ; eye(m) -Sb12])\full([SIh ; SIb]);
EHb = ([-Sh21 speye(m) ; speye(m) -Sb12])\([SIh ; SIb]);
%

%
[Ex,Ey,Ez] = deal(zeros(length(x),size(Ih,2)));
[Hx,Hy,Hz] = deal(zeros(length(x),size(Ih,2)));

%

[Em,Hm] = deal(EHb(1:m,:),EHb(1+m:end,:));
P = MatS{3};
Q = MatS{4};
InvP = inv(P);
InvQ = inv(Q);
% try
%     dP = decomposition(P,'lu');
%     dQ = decomposition(Q,'lu');
% catch
%     dP = P;
%     dQ = Q;
% end

Vp = MatS{6};
EpMuz = MatS{7};
if sqrt(size(EpMuz,1)) == m/2
    InvEpz = reshape(EpMuz(:,1),m/2,m/2);
    InvMuz = reshape(EpMuz(:,2),m/2,m/2);
    Epx = reshape(EpMuz(:,3),m/2,m/2);
    Epy = reshape(EpMuz(:,4),m/2,m/2);
else
    InvEpz = EpMuz(:,1)*eye(m/2);
    InvMuz = EpMuz(:,2)*eye(m/2);   
    Epx = EpMuz(:,3)*eye(m/2);   
    Epy = EpMuz(:,4)*eye(m/2);   
end
%
[MatSb1,MatSh1] = deal(MatS);

% Apodisation avec Hamming

Eps = .9;
TfApod = Apod(BetaX-BetaX0,BetaY-BetaY0,Eps);
%TfApod = ones(size(BetaX(:)));
%
MatQP = [Q , Q; P -P];
InvMatQP = 0.5*[InvQ , InvP; InvQ -InvP];
C = InvMatQP*EHb;
%
if issparse(InvMatQP), CondM = 1/condest(InvMatQP); else, CondM = rcond(InvMatQP); end

MatBetaX = repmat(reshape(BetaX(:),1,length(BetaX(:))),length(BetaX(:)),1);
MatBetaY = repmat(reshape(BetaY(:),1,length(BetaY(:))),length(BetaY(:)),1);
[MatPzBetaX,MatPzBetaY] = deal(InvEpz.*MatBetaX,InvEpz.*MatBetaY);
[MatQzBetaX,MatQzBetaY] = deal(InvMuz.*MatBetaX,InvMuz.*MatBetaY);

if isreal(MatPzBetaX), MatPzBetaX = complex(MatPzBetaX); end
if isreal(MatPzBetaY), MatPzBetaY = complex(MatPzBetaY); end
if isreal(MatQzBetaX), MatQzBetaX = complex(MatQzBetaX); end
if isreal(MatQzBetaY), MatQzBetaY = complex(MatQzBetaY); end

%[MatPzBetaX,MatPzBetaY] = deal(InvEpz*diag(BetaX),InvEpz*diag(BetaY));
%[MatQzBetaX,MatQzBetaY] = deal(InvMuz*diag(BetaX),InvMuz*diag(BetaY));

Phase = [];
epsz = max(abs(Mesh.CoorN(:)))/1e6;
xi = unique(Mesh.CoorV(:,1));
yi = unique(Mesh.CoorV(:,2));

while ~isempty(z)
    Pn = find(abs(z0 - min(z)) <= epsz);
    %Phase = exp(1i*Mesh.CoorN(Pn,1)*BetaX + 1i*Mesh.CoorN(Pn,2)*BetaY);
    %if isempty(Phase) || length(find(Pn)) ~= length(Phase) || (norm(x(Pn)-x0)>eps && norm(y(Pn)-y0)>eps)
        Phase = exp(1i*x(Pn)*BetaX + 1i*y(Pn)*BetaY);
    %end
    %
    %
    if max(Mesh.Nsd)>1 && length(unique(x(1:end-2))) ~= 1
        %
        for ki = 1:length(yi)
            Pe = find(abs(Mesh.CoorV(:,2)-yi(ki))<epsz);
            % 
            Pen = unique(Mesh.Cn(Pe,:));
            xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
            xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
            yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
            in = inpolygon(x(Pn),y(Pn),xd(:),yd(:));
            % test isempty(in)
            NumSd = unique(Mesh.Nsd(Pe));
            if length(unique(Data.Indice(NumSd))) == 1
                Ex(Pn(in),:) = Phase(in,:)*(TfApod.*Em(1:end/2,:));
            else
                Dx = Phase*(TfApod.*(Epx*Em(1:end/2,:)));
                for ke = 1:length(Pe)
                    ie = Pe(ke);
                    kd = Mesh.Nsd(ie);
                    Pen = unique(Mesh.Cn(ie,:));
                    xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
                    xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
                    yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
                    in = inpolygon(x(Pn),y(Pn),xd(:),yd(:));
                    %Pin = Pin & ~in;
                    Ex(Pn(in),:) = 1/Phys.CaractEps(1,1,kd)*Dx(in,:)/Phys.K0;
                end
            end
        end
    else
        Ex(Pn,:) = Phase*(TfApod.*Em(1:end/2,:));
    end
    %    %
    if max(Mesh.Nsd)>1 && length(unique(y(1:end-2))) ~= 1
        %
        for ki = 1:length(xi)
            Pe = find(abs(Mesh.CoorV(:,1)-xi(ki))<epsz);
            % 
            Pen = unique(Mesh.Cn(Pe,:));
            xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
            xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
            yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
            in = inpolygon(x(Pn),y(Pn),xd(:),yd(:));
            % test isempty(in)
            NumSd = unique(Mesh.Nsd(Pe));
            if length(unique(Data.Indice(NumSd))) == 1
                Ey(Pn(in),:) = Phase(in,:)*(TfApod.*Em(end/2+1:end,:));
            else
                Dy = Phase*(TfApod.*(Epy*Em(end/2+1:end,:)));
                for ke = 1:length(Pe)
                    ie = Pe(ke);
                    kd = Mesh.Nsd(ie);
                    Pen = unique(Mesh.Cn(ie,:));
                    xd = Mesh.CoorN(Pen,1); yd = Mesh.CoorN(Pen,2);
                    xd = [min(xd) max(xd) max(xd) min(xd) min(xd)];
                    yd = [min(yd) min(yd) max(yd) max(yd) min(yd)];
                    in = inpolygon(x(Pn),y(Pn),xd(:),yd(:));
                    %Pin = Pin & ~in;
                    Ey(Pn(in),:) = 1/Phys.CaractEps(2,2,kd)*Dy(in,:)/Phys.K0;
                end
            end
        end
    else
        Ey(Pn,:) = Phase*(TfApod.*Em(end/2+1:end,:));
    end
    %
    Hx(Pn,:) = Phase*(TfApod.*Hm(1:end/2,:));
    Hy(Pn,:) = Phase*(TfApod.*Hm(end/2+1:end,:));
    %Ezm = 1i*InvEpz*(diag(BetaX)*Hm(end/2+1:end) - diag(BetaY)*Hm(1:end/2));
    Ezm = 1i*(MatPzBetaX*Hm(end/2+1:end,:) - MatPzBetaY*Hm(1:end/2,:));
    Ez(Pn,:) = Phase*(TfApod.*Ezm);
    %Hzm = 1i*InvMuz*(diag(BetaX)*Em(end/2+1:end) - diag(BetaY)*Em(1:end/2));
    Hzm = 1i*(MatQzBetaX*Em(end/2+1:end,:) - MatQzBetaY*Em(1:end/2,:));
    Hz(Pn,:) = Phase*(TfApod.*Hzm);
    %
    x0 = x(Pn); y0 = y(Pn);
    %
    Pn = z > min(z);
    if isempty(z(Pn)), break; end
    hy = min(z(Pn))-zmin;
    %
    ExpVp = exp([Vp ;-Vp]*hy);
    %
    if 1/norm(ExpVp(end/2+1:end))*CondM < eps
        hy0 = hy + zmin - min(z0);

        InvCosH = repmat(reshape(1./cosh(Vp*hy0),1,length(Vp)),length(Vp),1);
        TanH = repmat(reshape(tanh(Vp*hy0),1,length(Vp)),length(Vp),1);
        MatSb1{1} = [(Q.*InvCosH)*InvQ , (Q.*TanH)*InvP; -(P.*TanH)*InvQ , (P.*InvCosH)*InvP];%[(Q.*InvCosH)/dQ , (Q.*TanH)/dP; -(P.*TanH)/dQ , (P.*InvCosH)/dP];
        %
        InvCosH = repmat(reshape(1./cosh(Vp*(hc-hy0)),1,length(Vp)),length(Vp),1);
        TanH = repmat(reshape(tanh(Vp*(hc-hy0)),1,length(Vp)),length(Vp),1);
        MatSh1{1} = [(Q.*InvCosH)*InvQ , (Q.*TanH)*InvP; -(P.*TanH)*InvQ , (P.*InvCosH)*InvP];%[(Q.*InvCosH)/dQ , (Q.*TanH)/dP; -(P.*TanH)/dQ , (P.*InvCosH)/dP];
        %
        %
        [MatSb1{2},MatSh1{2}] = deal(MatS{2});
        MatSh1 = ProdMatS(MatSh1,Sh);
        Sh1 = MatSh1{1};
        [Sh21,Sh22] = deal(Sh1(nh3+1:end,1:nh1),Sh1(nh3+1:end,nh1+1:end));
        MatSb1 = ProdMatS(Sb,MatSb1);
        Sb1 = MatSb1{1};
        [Sb11,Sb12] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end));
    
        SIh = Sh22*Ih(Pui,:);
        if size(SIh,2) == 0, SIh = zeros(size(SIh,1),size(Ih,2)); end
        SIb = Sb11*Ib(Pdi,:);
        if size(SIb,2) == 0, SIb = zeros(size(SIb,1),size(Ih,2)); end
    
        EHm = full([-Sh21 eye(m) ; eye(m) -Sb12])\full([SIh ; SIb]);
        C = InvMatQP*EHm;
        zmin = min(z(Pn));

    else
    %
        %ExpVp = exp([Vp ;-Vp]*hy);
        EHm = MatQP*(ExpVp.*C);
    end
    [Em,Hm] = deal(EHm(1:m,:),EHm(1+m:end,:));
    %
    z = z(Pn);
end

  
%norm(EHm-EHh)

% attention au clonage
if sum(Data.Sym) ~= 4
    VectE = [Ex(:,1) Ey(:,1) Ez(:,1)]*sqrt(Phys.Z0);
    VectH = [Hx(:,1) Hy(:,1) Hz(:,1)]/(1i*sqrt(Phys.Z0));
else
    VectE = [Ex(:,1) Ey(:,1) Ez(:,1) Ex(:,2) Ey(:,2) Ez(:,2)]*sqrt(Phys.Z0);
    VectH = [Hx(:,1) Hy(:,1) Hz(:,1) Hx(:,2) Hy(:,2) Hz(:,2)]/(1i*sqrt(Phys.Z0));
end
Vect = [VectE VectH];

end

%%
function TfApod = Apod(BetaX,BetaY,Eps)

%Eps = 0.75;

if length(BetaX(:))==1 && length(BetaY(:))==1, TfApod = 1; return; end

TfApodX = ones(size(BetaX(:)));
alphaX = BetaX*2*pi/(max(BetaX(:))-min(BetaX(:)));

Per = 2*Eps*pi; alpha0 = (1-Eps)*pi;
Px = find(alphaX>(1-Eps)*max(alphaX(:)));
%TfApodX(Px) = (0.5-0.5*cos(alphaX(Px)/Eps));
TfApodX(Px) = .5+.5*cos(2*pi/Per*(alphaX(Px)-alpha0));
Px = find(alphaX<-(1-Eps)*max(alphaX(:)));
%TfApodX(Px) = (0.5-0.5*cos(alphaX(Px)/Eps));
TfApodX(Px) = .5+.5*cos(2*pi/Per*(alphaX(Px)+alpha0));
%
TfApodY = ones(size(BetaY(:)));
alphaY = BetaY'*2*pi/(max(BetaY(:))-min(BetaY(:)));

Py = find(alphaY>(1-Eps)*max(alphaY(:)));
%TfApodY(Py) = (0.5-0.5*cos(alphaY(Py)/Eps));
TfApodY(Py) = .5+.5*cos(2*pi/Per*(alphaY(Py)-alpha0));
Py = find(alphaY<-(1-Eps)*max(alphaY(:)));
%TfApodY(Py) = (0.5-0.5*cos(alphaY(Py)/Eps));
TfApodY(Py) = .5+.5*cos(2*pi/Per*(alphaY(Py)+alpha0));


TfApod = TfApodX.*TfApodY;


end
