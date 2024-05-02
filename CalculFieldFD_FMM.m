function [E,H] = CalculFieldFD_FMM(Data,Mesh,Phys,Sb,MatS,Sh,X,Y,Z)

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

% Date of the latest version : 24 January 2023
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
            Mesh(k) = MoveMesh(utilMesh2(Mesh(k)),[0 0 max(Mesh(k-1).CoorN(:,end))]);
        else
            Mesh(k) = MoveMesh(utilMesh2(Mesh(k)),[0 max(Mesh(k-1).CoorN(:,end))]);
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
    error('For Field Calculation, symmetry is not considered');
end
%
% Numéro Couche calcul FD
if isfield(Data,'Nsub'), Nsub = cell2mat({Data.Nsub}); k0 = find(Nsub~=0); end
%
[E,H] = deal(nan(length(x(:)),3*size(MatS{min(k0),6},2)/2));

a0 = max(Mesh(1).CoorN(:,1))-min(Mesh(1).CoorN(:,1));
b0 = max(Mesh(1).CoorN(:,2))-min(Mesh(1).CoorN(:,2));



% Milieu bas
if min(z(:))<min(Mesh(1).CoorN(:,end))
    %
    h0 = min(Mesh(1).CoorN(:,end))-min(z(:));
    %
    if size(Mesh(1).CoorN,2) == 3
        Mesh0 = MeshLayer(a0,[],b0,[],h0,2,2,2);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 0 min(z(:))+min(Mesh(1).CoorN(:,end))]);
    elseif size(Mesh(1).CoorN,2) == 2
        Mesh0 = MeshLayer(a0,[],h0,2,2);
        Mesh0 = MoveMesh(utilMesh2(Mesh0),[0 min(z(:))+min(Mesh(1).CoorN(:,end))]);
    end
    %
    Data0 = SetData(Data(1),'Indice',Data(1).nb,'Nsub',0);
    Phys0 = CaractMat(Mesh0,Data0);
    S0 = CalculMatS(Data0,Mesh0,Phys0); % milieu bas 
    if isfield(Mesh(1),'xv'), Mesh0.xv = []; Mesh0.yv = []; end
    Mesh = [Mesh0 Mesh];
    MatS = [S0; MatS];
    Phys = [Phys0 Phys];
    Data = [Data0 Data];
    k0 = k0+1;
    [r,t,CoefD,R,T,Em,Hm] = CalculCoefRT(Sb,MatS,Sh,k0);
    MatS{k0,6} = [Em Hm];
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
    Data0 = SetData(Data(1),'Indice',Data(1).nh,'Nsub',0);
    Phys0 = CaractMat(Mesh0,Data0);
    S0 = CalculMatS(Data0,Mesh0,Phys0); % milieu haut 
    if isfield(Mesh(1),'xv'), Mesh0.xv = []; Mesh0.yv = []; end
    Mesh = [Mesh Mesh0];
    MatS = [MatS; S0];
    Phys = [Phys Phys0];

    Data = [Data Data0];
end
%
% Si on ajoute des couches il faut recalculer Em Hm du réseau
%
%tic
% Numéro Couche calcul FD
%if isfield(Data,'Nsub'), Nsub = cell2mat({Data.Nsub}); k0 = find(Nsub~=0); end


for k = 1:length(Mesh)
    % Grille de calcul
    Pn = find(z(:)>=min(Mesh(k).CoorN(:,end)) & z(:)<=max(Mesh(k).CoorN(:,end)));
    %
    if ~isempty(Pn)
        x0 = x(Pn); y0 = y(Pn); z0 = z(Pn);
        %
        if all(k~= k0)
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
        end
        %
        % Mise à jour de Sb1 et Sh1
        % Calcul du champ par cellule
        if any(k == k0)
            [VectE,VectH] = FieldFD_MMF3D(Data(k),Mesh(k),Phys(k),MatS(k,:),[x0(:) y0(:) z0(:)]);
        elseif k < min(k0)
            if k == 1, Sb1 = Sb; else, Sb1 = ProdMatS(Sb,ProdMatS(MatS(1:k-1,:))); end
            if k+1 == min(k0) , Sh1 = []; else, Sh1 = ProdMatS(MatS(k+1:min(k0)-1,:)); end
            [VectE,VectH] = FieldFD_MMF_b(Data(k),Mesh(k),Phys(k),Sb1,MatS(k,:),Sh1,[x0(:) y0(:) z0(:)],MatS{min(k0),6});
    
        elseif k > max(k0)
            if k == length(Mesh), Sh1 = Sh; else, Sh1 = ProdMatS(ProdMatS(MatS(k+1:end,:)),Sh); end
            if k-1 == max(k0), Sb1 = []; else, Sb1 = ProdMatS(MatS(max(k0)+1:k-1,:)); end
            [VectE,VectH] = FieldFD_MMF_h(Data(k),Mesh(k),Phys(k),Sb1,MatS(k,:),Sh1,[x0(:) y0(:) z0(:)],MatS{max(k0),6});
        end
        % 
        E(Pn,1:end) = VectE(1:length(Pn),:);
        H(Pn,1:end) = VectH(1:length(Pn),:);
    end
end
%toc

end

% ------------------------------------------------------------------------%




%% -------------------------------------------------------------------------
% Structure 3D
% -------------------------------------------------------------------------

function [VectE,VectH,Vect] = FieldFD_MMF3D(Data,Mesh,Phys,MatS,TabCoor)

if isempty(TabCoor)
    [x,y,z] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2),Mesh.CoorN(:,3));
else
    Pn = TabCoor(:,3)>=min(Mesh.CoorN(:,3)) & TabCoor(:,3)<=max(Mesh.CoorN(:,3));
    [x,y,z] = deal(TabCoor(Pn,1),TabCoor(Pn,2),TabCoor(Pn,3));
end
z0 = z;
%
Np = Data.Nsub;

z1 = linspace(min(Mesh.CoorN(:,3)),max(Mesh.CoorN(:,3)),Np+1);

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
if iscell(MatS{1}), MatS{1} = cell2mat(MatS{1}); end
m = size(MatS{1},1); 
%
%
[Ex,Ey,Ez] = deal(zeros(length(x),2));
[Hx,Hy,Hz] = deal(zeros(length(x),2));
%

EpMuz = MatS{7};
if sqrt(size(EpMuz,1)) == m/2
    InvEpz = reshape(EpMuz(:,1),m/2,m/2);
    InvMuz = reshape(EpMuz(:,2),m/2,m/2);
else
    InvEpz = EpMuz(:,1)*eye(m/2);
    InvMuz = EpMuz(:,2)*eye(m/2);   
end
%

%    TfApod = ones(size(BetaX(:)));
Eps = 0.25;
TfApod = Apod(BetaX,BetaY,Eps);
%
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
%Phase = exp(1i*x(:)*BetaX + 1i*y(:)*BetaY);
Phase = [];

while ~isempty(z)
    Pn = abs(z0 - min(z)) <= max(abs(Mesh.CoorN(:)))/1e6;
    Pk = find((z1-min(z))<=0);
    k = Pk(end);
    [Em1,Hm1] = deal(MatS{6}((1:m)+(k-1)*m,1:2),MatS{6}((1:m)+(k-1)*m,3:4)); % mettre à jour Em Hm 

    if k == length(z1)
        Alpha = 0; 
        [Em2,Hm2] = deal(0,0);
    else 
        Alpha = (min(z)-z1(k))/(z1(k+1)-z1(k)); 
        [Em2,Hm2] = deal(MatS{6}((1:m)+k*m,1:2),MatS{6}((1:m)+k*m,3:4)); % mettre à jour Em Hm 
    end
    %
    Em = (1-Alpha)*Em1 + Alpha*Em2;
    Hm = (1-Alpha)*Hm1 + Alpha*Hm2;
    
    if isempty(Phase) || (norm(x(Pn)-x0)>eps && norm(y(Pn)-y0)>eps)
        Phase = exp(1i*x(Pn)*BetaX + 1i*y(Pn)*BetaY);
    end
    %
    Ex(Pn,:) = Phase*(TfApod.*Em(1:end/2,:));
    Ey(Pn,:) = Phase*(TfApod.*Em(end/2+1:end,:));
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
    Pn = z > min(z);
    z = z(Pn);
end

% attention au clonage
VectE = [Ex(:,1) Ey(:,1) Ez(:,1) Ex(:,2) Ey(:,2) Ez(:,2)]*sqrt(Phys.Z0);
VectH = [Hx(:,1) Hy(:,1) Hz(:,1) Hx(:,2) Hy(:,2) Hz(:,2)]/(1i*sqrt(Phys.Z0));
Vect = [VectE VectH];

end

%% -------------------------------------------------------------------------
% Structure 3D milieu bas
% -------------------------------------------------------------------------

function [VectE,VectH,Vect] = FieldFD_MMF_b(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor,EH)

if isempty(TabCoor)
    [x,y,z] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2),Mesh.CoorN(:,3));
    z0 = z;
else
    Pn = TabCoor(:,3)>=min(Mesh.CoorN(:,3)) & TabCoor(:,3)<=max(Mesh.CoorN(:,3));
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
Ib = zeros(2*m,2); % % dim 2 A vérifier

mx = m/2;
%

Ib((mx-1)/2+1+mx,1) = 1;
Ib((mx-1)/2+1,2) = 1;
%
[Pdi,Pd] = deal(Sb{7},Sb{8});

Sb1 = Sb{1};
[nb1,nb2,nb3,nb4] = deal(length(Pdi),m,m,length(Pd));
[Sb11,Sb12] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end));
%
[Eb,Hb] = deal(EH(1:m,1:2),EH(1:m,3:4));

if isempty(Sh), Sh1 = MatS{1}; else, MatSh1 = ProdMatS(MatS,Sh); Sh1 = MatSh1{1}; end
[Sh21,Sh22] = deal(Sh1(1+m:end,1:m),Sh1(1+m:end,1+m:end));
SIh = Sh22*Hb;

SIb = Sb11*Ib(Pdi,:);
if size(SIb,2) == 0, SIb = zeros(size(SIb,1),2); end


EHb = full([-Sh21 eye(m) ; eye(m) -Sb12])\full([SIh ; SIb]);
%
%
%

%
[Ex,Ey,Ez] = deal(zeros(length(x),2));
[Hx,Hy,Hz] = deal(zeros(length(x),2));

%

[Em,Hm] = deal(EHb(1:m,:),EHb(1+m:end,:));
P = MatS{3};
Q = MatS{4};
InvP = inv(P);
InvQ = inv(Q);
%
Vp = MatS{6};
EpMuz = MatS{7};
if sqrt(size(EpMuz,1)) == m/2
    InvEpz = reshape(EpMuz(:,1),m/2,m/2);
    InvMuz = reshape(EpMuz(:,2),m/2,m/2);
else
    InvEpz = EpMuz(:,1)*eye(m/2);
    InvMuz = EpMuz(:,2)*eye(m/2);   
end
%
[MatSb1,MatSh1] = deal(MatS);

% Apodisation avec Hamming

Eps = 0.25;
TfApod = Apod(BetaX,BetaY,Eps);
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

while ~isempty(z)
    Pn = abs(z0 - min(z)) <= max(abs(Mesh.CoorN(:)))/1e6;
    %Phase = exp(1i*Mesh.CoorN(Pn,1)*BetaX + 1i*Mesh.CoorN(Pn,2)*BetaY);
    if isempty(Phase) || (norm(x(Pn)-x0)>eps && norm(y(Pn)-y0)>eps)
        Phase = exp(1i*x(Pn)*BetaX + 1i*y(Pn)*BetaY);
    end
    %
    Ex(Pn,:) = Phase*(TfApod.*Em(1:end/2,:));
    Ey(Pn,:) = Phase*(TfApod.*Em(end/2+1:end,:));
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
        %MatSh1 = ProdMatS(MatSh1,Sh);
        %Sh1 = MatSh1{1};
        if isempty(Sh), Sh1 = MatSh1{1}; else, MatSh1 = ProdMatS(MatSh1,Sh); Sh1 = MatSh1{1}; end
        [Sh21,Sh22] = deal(Sh1(1+m:end,1:m),Sh1(1+m:end,1+m:end));
        SIh = Sh22*Hb;


        MatSb1 = ProdMatS(Sb,MatSb1);
        Sb1 = MatSb1{1};
        [Sb11,Sb12] = deal(Sb1(1:nb3,1:nb1),Sb1(1:nb3,nb1+1:end));
    
        %SIh = Sh22*Ih(Pui);
        %if size(SIh,2) == 0, SIh = zeros(size(SIh,1),1); end
        SIb = Sb11*Ib(Pdi,:);
        if size(SIb,2) == 0, SIb = zeros(size(SIb,1),2); end
    
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
VectE = [Ex(:,1) Ey(:,1) Ez(:,1) Ex(:,2) Ey(:,2) Ez(:,2)]*sqrt(Phys.Z0);
VectH = [Hx(:,1) Hy(:,1) Hz(:,1) Hx(:,2) Hy(:,2) Hz(:,2)]/(1i*sqrt(Phys.Z0));
Vect = [VectE VectH];

end

%% -------------------------------------------------------------------------
% Structure 3D milieu haut
% -------------------------------------------------------------------------

function [VectE,VectH,Vect] = FieldFD_MMF_h(Data,Mesh,Phys,Sb,MatS,Sh,TabCoor,EH)

if isempty(TabCoor)
    [x,y,z] = deal(Mesh.CoorN(:,1),Mesh.CoorN(:,2),Mesh.CoorN(:,3));
    z0 = z;
else
    Pn = TabCoor(:,3)>=min(Mesh.CoorN(:,3)) & TabCoor(:,3)<=max(Mesh.CoorN(:,3));
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
m = size(Sh{3},1); 
Ih = zeros(2*m,2); % % dim 2 A vérifier

mx = m/2;

%
Ih((mx-1)/2+1+3*mx,1) = 1;
Ih((mx-1)/2+1+2*mx,2) = 1;

%
[Pui,Pu] = deal(Sh{7},Sh{8});

MatSh1 = ProdMatS(MatS,Sh);
Sh1 = MatSh1{1};
[nh1,nh2,nh3,nh4] = deal(m,length(Pui),length(Pu),m);
[Sh21,Sh22] = deal(Sh1(nh3+1:end,1:nh1),Sh1(nh3+1:end,nh1+1:end));
SIh = Sh22*Ih(Pui,:);  % A vérifier
if size(SIh,2) == 0, SIh = zeros(size(SIh,1),2); end
%
[Eh,Hh] = deal(EH(end-m+1:end,1:2),EH(end-m+1:end,3:4));

if isempty(Sb)
    EHb = [Eh;Hh];
else
    Sb1 = Sb{1};
    [Sb11,Sb12] = deal(Sb1(1:m,1:m),Sb1(1:m,m+1:end));
    SIb = Sb11*Eh;
    EHb = full([-Sh21 eye(m) ; eye(m) -Sb12])\full([SIh ; SIb]);
end
%
%
[Ex,Ey,Ez] = deal(zeros(length(x),2));
[Hx,Hy,Hz] = deal(zeros(length(x),2));

%

[Em,Hm] = deal(EHb(1:m,:),EHb(1+m:end,:));
P = MatS{3};
Q = MatS{4};
InvP = inv(P);
InvQ = inv(Q);
%

Vp = MatS{6};
EpMuz = MatS{7};
if sqrt(size(EpMuz,1)) == m/2
    InvEpz = reshape(EpMuz(:,1),m/2,m/2);
    InvMuz = reshape(EpMuz(:,2),m/2,m/2);
else
    InvEpz = EpMuz(:,1)*eye(m/2);
    InvMuz = EpMuz(:,2)*eye(m/2);   
end
%
[MatSb1,MatSh1] = deal(MatS);

% Apodisation avec Hamming
Eps = 0.25;
TfApod = Apod(BetaX,BetaY,Eps);

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

while ~isempty(z)
    Pn = abs(z0 - min(z)) <= max(abs(Mesh.CoorN(:)))/1e6;
    %Phase = exp(1i*Mesh.CoorN(Pn,1)*BetaX + 1i*Mesh.CoorN(Pn,2)*BetaY);
    if isempty(Phase) || (norm(x(Pn)-x0)>eps && norm(y(Pn)-y0)>eps)
        Phase = exp(1i*x(Pn)*BetaX + 1i*y(Pn)*BetaY);
    end
    %
    Ex(Pn,:) = Phase*(TfApod.*Em(1:end/2,:));
    Ey(Pn,:) = Phase*(TfApod.*Em(end/2+1:end,:));
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
        SIh = Sh22*Ih(Pui,:);
        if size(SIh,2) == 0, SIh = zeros(size(SIh,1),2); end

        if ~isempty(Sb), MatSb1 = ProdMatS(Sb,MatSb1); end
        Sb1 = MatSb1{1};
        [Sb11,Sb12] = deal(Sb1(1:m,1:m),Sb1(1:m,m+1:end));
    
        SIb = Sb11*Eh;
    
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
VectE = [Ex(:,1) Ey(:,1) Ez(:,1) Ex(:,2) Ey(:,2) Ez(:,2)]*sqrt(Phys.Z0);
VectH = [Hx(:,1) Hy(:,1) Hz(:,1) Hx(:,2) Hy(:,2) Hz(:,2)]/(1i*sqrt(Phys.Z0));
Vect = [VectE VectH];

end

%%
function TfApod = Apod(BetaX,BetaY,Eps)

%Eps = 0.25;

TfApodX = ones(size(BetaX(:)));
alphaX = BetaX*2*pi/(max(BetaX(:))-min(BetaX(:)));

Px = find(alphaX>(1-Eps)*max(alphaX(:)));
TfApodX(Px) = (0.5-0.5*cos(alphaX(Px)/Eps));
Px = find(alphaX<-(1-Eps)*max(alphaX(:)));
TfApodX(Px) = (0.5-0.5*cos(alphaX(Px)/Eps));
%
TfApodY = ones(size(BetaY(:)));
alphaY = BetaY'*2*pi/(max(BetaY(:))-min(BetaY(:)));

Py = find(alphaY>(1-Eps)*max(alphaY(:)));
TfApodY(Py) = (0.5-0.5*cos(alphaY(Py)/Eps));
Py = find(alphaY<-(1-Eps)*max(alphaY(:)));
TfApodY(Py) = (0.5-0.5*cos(alphaY(Py)/Eps));

TfApod = TfApodX.*TfApodY;

end

