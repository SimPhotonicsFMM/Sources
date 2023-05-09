function Phys = CaractMat(Mesh,Data,Kx,Ky)

% CaractMat 
%   Description ofmaterial properperties 
%   This function is also used for FEM problem 
%
% Syntax
%   Phys = CaractMat();  % c0, Eps0, Mu0, Z0
%   Phys = CaractMat(Mesh,Data);
%   Phys = CaractMat(Mesh,Data,Kx); 
%   Phys = CaractMat(Mesh,Data,Kx,Ky);
%
% Description
%   Data : Data of the priblem
%   Mesh : Mesh of the structure (CoorN, Cn, Nsd,...CoorA, Ca, ExtAr, )
%
%   Phys : Material properties (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%

% Date of the latest version : 06 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if nargin == 0
    Phys.c0 = 2.99792458e+14;        % µm/s
    Phys.Mu0 = 4*pi*1e-7;            % Perméabilité du vide
    Phys.Eps0 = 1/(36*pi*1e9);       % Permittivité du vide
    Phys.Z0 = sqrt(Phys.Mu0/Phys.Eps0);
    return
end

%
switch size(Mesh(1).CoorN,2)
    case 1
        Phys = CaractMat2D(Mesh,Data);
    case 2
        if isfield(Data(1),'DofType') && strcmp(Data(1).DofType,'Edge')
            for kc = 1:length(Mesh) 
                Phys(kc) = CaractMat3D(Mesh(kc),Data(kc));
            end

        else
            for kc = 1:length(Mesh) 
                Phys(kc) = CaractMat2D(Mesh(kc),Data(kc));
            end
        end
    case 3
           for kc = 1:length(Mesh) 
                Phys(kc) = CaractMat3D(Mesh(kc),Data(kc));
           end
           %Phys.TabIndice = zeros(size(Mesh.Nsd));
        %for k = 1:max(Mesh.Nsd), Phys.TabIndice(Mesh.Nsd == k) = Data.Indice(k); end
        

end
% Mise à jour de Kx et Ky
for kc = 1:length(Mesh) 
    if nargin >= 3, Phys(kc).Kx = Kx; end
    if nargin >= 4, Phys(kc).Ky = Ky; end
end

if isfield(Data,'Sigma')
    for k = 1:length(Data.Sigma)
        Phys.Sigma(1:size(Mesh.Cn,1),k) = Data.Sigma(k); 
    end
end

% Caractéristiques thermiques et mécaniques

if ~isfield(Data,'Fv'), return, end


Phys.TabFs = Data.TabFs;
Phys.Fv = Data.Fv;

% En coordonnées sphériques 1D ou cylindriques 2D
NDim = str2double(Data.Dim(1));
Phys.Cond = zeros(size(Mesh.Cn,1),NDim);

for k = 1:max(Mesh.Nsd)
    if strcmp(Data.Dim,'1DSpherical')
        Rayon = Mesh.CoorA;
        Phys.Cond(Mesh.Nsd == k) = Data.Cond(k)*4*pi*Rayon(Mesh.Nsd == k).^2;
        Phys.Cp(Mesh.Nsd == k) = Data.Cp(k)*4*pi*Rayon(Mesh.Nsd == k).^2;
        Phys.Qv(Mesh.Nsd == k) = Data.Qv(k)*4*pi*Rayon(Mesh.Nsd == k).^2;
    elseif strcmp(Data.Dim,'1DAxi') 
        Rayon = Mesh.CoorA(:,1);
        Phys.Cond(Mesh.Nsd == k) = Data.Cond(k)*2*pi*Rayon(Mesh.Nsd == k);
        Phys.Cp(Mesh.Nsd == k) = Data.Cp(k)*2*pi*Rayon(Mesh.Nsd == k);
        Phys.Qv(Mesh.Nsd == k) = Data.Qv(k)*2*pi*Rayon(Mesh.Nsd == k);
    elseif strcmp(Data.Dim,'2DAxi')
        Rayon = Mesh.CoorF(:,1);
        Phys.Cond(Mesh.Nsd == k) = Data.Cond(k)*2*pi*Rayon(Mesh.Nsd == k);
        Phys.Cp(Mesh.Nsd == k) = Data.Cp(k)*2*pi*Rayon(Mesh.Nsd == k);
        Phys.Qv(Mesh.Nsd == k) = Data.Qv(k)*2*pi*Rayon(Mesh.Nsd == k);
    else
        
        if numel(Data.Cond) == max(Mesh.Nsd)
            Phys.Cond(Mesh.Nsd == k,1:NDim) = Data.Cond(k);
        else
            for kk = 1:NDim, Phys.Cond(Mesh.Nsd == k,kk) = Data.Cond(k,kk); end
        end
        Phys.Cp(Mesh.Nsd == k) = Data.Cp(k);
        Phys.Qv(Mesh.Nsd == k) = Data.Qv(k);
        Phys.ModYoung(Mesh.Nsd == k) = Data.ModYoung(k);
        Phys.CoefPoisson(Mesh.Nsd == k) = Data.CoefPoisson(k);
        Phys.MasseVol(Mesh.Nsd == k) = Data.MasseVol(k);
        Phys.Alpha(Mesh.Nsd == k) = Data.Alpha(k);
    end  
end

% Phys.TabHc = Data.TabHc;
% Phys.TabHr = Data.TabHr*5.6703e-8;
% Phys.TabQs = Data.TabQs;


if strcmp(Data.Dim,'1DSpherical')
    Rayon = Mesh.CoorN;
    Phys.TabHc(1) = Data.TabHc(1)*4*pi*min(Rayon)^2;
    Phys.TabHc(2) = Data.TabHc(2)*4*pi*max(Rayon)^2;
    %
    Phys.TabHr(1) = Data.TabHr(1)*4*pi*min(Rayon)^2*5.6703e-8;
    Phys.TabHr(2) = Data.TabHr(2)*4*pi*max(Rayon)^2*5.6703e-8;
    %
    Phys.TabQs(1) = Data.TabQs(1)*4*pi*min(Rayon)^2;
    Phys.TabQs(2) = Data.TabQs(2)*4*pi*max(Rayon)^2;
elseif strcmp(Data.Dim,'1DAxi')
    Rayon = Mesh.CoorN;
    Phys.TabHc(1) = Data.TabHc(1)*2*pi*min(Rayon);
    Phys.TabHc(2) = Data.TabHc(2)*2*pi*max(Rayon);
    %
    Phys.TabHr(1) = Data.TabHr(1)*2*pi*min(Rayon)*5.6703e-8;
    Phys.TabHr(2) = Data.TabHr(2)*2*pi*min(Rayon)*5.6703e-8;
    %
    Phys.TabQs(1) = Data.TabQs(1)*2*pi*min(Rayon);
    Phys.TabQs(2) = Data.TabQs(2)*2*pi*min(Rayon);
elseif strcmp(Data.Dim,'1D')
    Phys.TabHc = Data.TabHc;
    Phys.TabHr = Data.TabHr*5.6703e-8;
    Phys.TabQs = Data.TabQs;
elseif strcmp(Data.Dim,'2DAxi')
    Rayon = Mesh.CoorA(:,1);
    RayonA = sqrt(sum(Mesh.CoorA.^2,2));
    [Phys.TabHc,Phys.TabHr,Phys.TabQs] = deal(zeros(length(Mesh.CoorA),1));
    % Rmin (ligne)
    Ps = abs(Rayon - min(Rayon))<= max(abs(Mesh.CoorA(:,1)))/1e4;
    Phys.TabHc(Ps) = Data.TabHc(1)*2*pi*Rayon(Ps);
    Phys.TabHr(Ps) = Data.TabHr(1)*2*pi*Rayon(Ps)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1)*2*pi*Rayon(Ps);
    % Rmax (Ligne)
    Ps = abs(Rayon - max(Rayon))<= max(abs(Mesh.CoorA(:,1)))/1e4;
    Phys.TabHc(Ps) = Data.TabHc(1+size(Data.TabHc,1))*2*pi*Rayon(Ps);
    Phys.TabHr(Ps) = Data.TabHr(1+size(Data.TabHc,1))*2*pi*Rayon(Ps)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1+size(Data.TabHc,1))*2*pi*Rayon(Ps);
    % Zmin
    Ps = abs(Mesh.CoorA(:,2) - min(Mesh.CoorA(:,2)))<= max(abs(Mesh.CoorA(:,1)))/1e4;
    Phys.TabHc(Ps) = Data.TabHc(2)*2*pi*Rayon(Ps);
    Phys.TabHr(Ps) = Data.TabHr(2)*2*pi*Rayon(Ps)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2)*2*pi*Rayon(Ps);
    % Zmax
    Ps = abs(Mesh.CoorA(:,2) - max(Mesh.CoorA(:,2)))<= max(abs(Mesh.CoorA(:,1)))/1e4;
    Phys.TabHc(Ps) = Data.TabHc(2+size(Data.TabHc,1))*2*pi*Rayon(Ps);
    Phys.TabHr(Ps) = Data.TabHr(2+size(Data.TabHc,1))*2*pi*Rayon(Ps)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2+size(Data.TabHc,1))*2*pi*Rayon(Ps);
    % Contour circulaire
    if size(Data.TabHc,1) == 3,
        % Rmin (cercle)
        Ps = abs(RayonA - min(RayonA))<= max(abs(Mesh.CoorA(:,1)))/1e4;
        Phys.TabHc(Ps) = Data.TabHc(3)*2*pi*Rayon(Ps);
        Phys.TabHr(Ps) = Data.TabHr(3)*2*pi*Rayon(Ps)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(3)*2*pi*Rayon(Ps);
        % Rmax (cercle)
        %Ps = abs(RayonA - max(RayonA))<= max(abs(Mesh.CoorA(:,1)))/1e4;
        Pa = zeros(length(Mesh.CoorA),1);
        for ka = 1:length(Mesh.CoorA), Pa(ka) = sum(sum(ismember(abs(Mesh.Ca),ka),2)); end
        Ps = find(Pa==1 & Mesh.CoorA(:,1)>0); % à comparer rayon ?
        Phys.TabHc(Ps) = Data.TabHc(6)*2*pi*Rayon(Ps);
        Phys.TabHr(Ps) = Data.TabHr(6)*2*pi*Rayon(Ps)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(6)*2*pi*Rayon(Ps);
    end
elseif strcmp(Data.Dim,'2D')
    RayonA = sqrt(sum(Mesh.CoorA.^2,2));
    [Phys.TabHc,Phys.TabHr,Phys.TabQs] = deal(zeros(length(Mesh.CoorA),1));
    % Rmin (cercle)
    if size(Data.TabHc,1) == 3,
        Ps = abs(RayonA - min(RayonA))<= max(abs(Mesh.CoorA(:,1)))/1e4;
        Phys.TabHc(Ps) = Data.TabHc(3);
        Phys.TabHr(Ps) = Data.TabHr(3)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(3);
    % Rmax (cercle)
        Ps = abs(RayonA - max(RayonA))<= max(abs(Mesh.CoorA(:,1)))/1e4;
        Phys.TabHc(Ps) = Data.TabHc(6);
        Phys.TabHr(Ps) = Data.TabHr(6)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(6);
    end
    % Xmin 
    Ps = abs(Mesh.CoorA(:,1) - min(Mesh.CoorA(:,1)))<= max(abs(Mesh.CoorA(:,1)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(1);
    Phys.TabHr(Ps) = Data.TabHr(1)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1);
    % Xmax 
    Ps = abs(Mesh.CoorA(:,1) - max(Mesh.CoorA(:,1)))<= max(abs(Mesh.CoorA(:,1)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(1+size(Data.TabHc,1));
    Phys.TabHr(Ps) = Data.TabHr(1+size(Data.TabHc,1))*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1+size(Data.TabHc,1));
    % Ymin
    Ps = abs(Mesh.CoorA(:,2) - min(Mesh.CoorA(:,2)))<= max(abs(Mesh.CoorA(:,2)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(2);
    Phys.TabHr(Ps) = Data.TabHr(2)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2);
    % Ymax
    Ps = abs(Mesh.CoorA(:,2) - max(Mesh.CoorA(:,2)))<= max(abs(Mesh.CoorA(:,2)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(2+size(Data.TabHc,1));
    Phys.TabHr(Ps) = Data.TabHr(2+size(Data.TabHc,1))*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2+size(Data.TabHc,1));
elseif strcmp(Data.Dim,'3D')
    RayonF = sqrt(sum(Mesh.CoorF.^2,2));
    [Phys.TabHc,Phys.TabHr,Phys.TabQs,Phys.TabFs] = deal(zeros(length(Mesh.CoorF),1));
    % Rmin (cercle)
    % Xmin
    Ps = abs(Mesh.CoorF(:,1) - min(Mesh.CoorF(:,1)))<= max(abs(Mesh.CoorA(:,1)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(1);
    Phys.TabHr(Ps) = Data.TabHr(1)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1);
    Phys.TabFs(Ps) = Data.TabFs(1);
    % Xmax 
    Ps = abs(Mesh.CoorF(:,1) - max(Mesh.CoorF(:,1)))<= max(abs(Mesh.CoorA(:,1)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(1+size(Data.TabHc,1));
    Phys.TabHr(Ps) = Data.TabHr(1+size(Data.TabHc,1))*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(1+size(Data.TabHc,1));
    Phys.TabFs(Ps) = Data.TabFs(1+size(Data.TabFs,1));
    % Ymin
    Ps = abs(Mesh.CoorF(:,2) - min(Mesh.CoorF(:,2)))<= max(abs(Mesh.CoorA(:,2)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(2);
    Phys.TabHr(Ps) = Data.TabHr(2)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2);
    Phys.TabFs(Ps) = Data.TabFs(2);
    % Ymax
    Ps = abs(Mesh.CoorF(:,2) - max(Mesh.CoorF(:,2)))<= max(abs(Mesh.CoorA(:,2)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(2+size(Data.TabHc,1));
    Phys.TabHr(Ps) = Data.TabHr(2+size(Data.TabHc,1))*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(2+size(Data.TabHc,1));
    Phys.TabFs(Ps) = Data.TabFs(2+size(Data.TabFs,1));
    % Zmin
    Ps = abs(Mesh.CoorF(:,3) - min(Mesh.CoorF(:,3)))<= max(abs(Mesh.CoorA(:,3)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(3);
    Phys.TabHr(Ps) = Data.TabHr(3)*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(3);
    Phys.TabFs(Ps) = Data.TabFs(3);
    % Zmax
    Ps = abs(Mesh.CoorF(:,3) - max(Mesh.CoorF(:,3)))<= max(abs(Mesh.CoorA(:,3)))/1e5;
    Phys.TabHc(Ps) = Data.TabHc(3+size(Data.TabHc,1));
    Phys.TabHr(Ps) = Data.TabHr(3+size(Data.TabHc,1))*5.6703e-8;
    Phys.TabQs(Ps) = Data.TabQs(3+size(Data.TabHc,1));
    Phys.TabFs(Ps) = Data.TabFs(3+size(Data.TabFs,1));
    if size(Data.TabHc,1) == 4
        Ps = abs(RayonF -  min(RayonF))<= max(abs(RayonF))/1e3;
        Phys.TabHc(Ps) = Data.TabHc(4);
        Phys.TabHr(Ps) = Data.TabHr(4)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(4);
    % Rmax (cercle)
        Ps = abs(RayonF - max(RayonF))<= max(abs(RayonF))/1e3;
        Phys.TabHc(Ps) = Data.TabHc(8);
        Phys.TabHr(Ps) = Data.TabHr(8)*5.6703e-8;
        Phys.TabQs(Ps) = Data.TabQs(8);
    end

end


end

%%

function Phys = CaractMat2D(Mesh,Data)

% CaractMat2D 
%   Description des matériaux et prise en compte des PML 
%
% Syntaxe
%   Phys = CaractMat2D(Mesh,Data);
%
% Description
%   Mesh : structure données du maillage (CoorN, Cn, CoorA, Ca, ...,Nsd )
%   Data : données du problème
%
%   Phys : structure caractéristiques des matériaux
%       Indice : indice de réfraction des matériaux
%       PmlX : Numéros des sous-domaines de PML en x
%       PmlY : Numéros des sous-domaines de PML en y
%
%   Date de la dernière version : 07 Décembre 2009
%   Auteur : Mondher Besbes (LCF / CNRS / IOGS)

Phys.PmlX = [];
Phys.PmlY = [];

if size(Mesh.CoorN,2) == 2
    if numel(Data.a) == 1
        if isfinite(Data.a), Phys.PmlX = unique(Mesh.Nsd(abs(Mesh.CoorF(:,1)) > Data.a)); end
    else
        Phys.PmlX = unique(Mesh.Nsd(Mesh.CoorF(:,1) < Data.a(1) | Mesh.CoorF(:,1) > Data.a(2)));
    end
    %
    if numel(Data.b) == 1
        if isfinite(Data.b), Phys.PmlY = unique(Mesh.Nsd(abs(Mesh.CoorF(:,2)) > Data.b)); end
    else
        Phys.PmlY = unique(Mesh.Nsd(Mesh.CoorF(:,2) < Data.b(1) | Mesh.CoorF(:,2) > Data.b(2)));
    end

elseif size(Mesh.CoorN,2) == 1
    if numel(Data.a) == 1
        if isfinite(Data.a), Phys.PmlX = unique(Mesh.Nsd(abs(Mesh.CoorA(:,1)) > Data.a)); end
    else
        Phys.PmlX = unique(Mesh.Nsd(Mesh.CoorA(:,1) < Data.a(1) | Mesh.CoorA(:,1) > Data.a(2)));
    end
end



Phys.TypePol = Data.TypePol;
Phys.CoefPml = Data.CoefPml;
Phys.Lambda0 = Data.Lambda0;
Phys.c0 = 2.99792458e+14;        % µm/s
Phys.Mu0 = 4*pi*1e-7;            % Perméabilité du vide
Phys.Eps0 = 1/(36*pi*1e9);       % Permittivité du vide
Phys.Z0 = sqrt(Phys.Mu0/Phys.Eps0);
Phys.K0 = 2*pi/Phys.Lambda0;
Phys.Indice = Data.Indice;
Phys.Sym = Data.Sym;

%
CaractNu = zeros(3,3,max(Mesh.Nsd));
CaractEps = zeros(3,3,max(Mesh.Nsd));
CaractNu1 = zeros(3,3,max(Mesh.Nsd));

switch Phys.TypePol
    case 0               % Polarisation TE
        for ic = 1:3, CaractNu(ic,ic,:) = 1; end
        for ic = 1:3, CaractNu1(ic,ic,:) = 1; end

        for ic = 1:3
            if iscell(Data.Indice)
                CaractEps(ic,ic,:) = Data.Indice{ic}.^2; 
            else
                CaractEps(ic,ic,:) = Data.Indice.^2; 
            end
        end
        MatPml = CaractNu;
    case 2               % Polarisation TM
        for ic = 1:3 
            CaractEps(ic,ic,:) = 1; 
                if iscell(Data.Indice)
                    CaractNu(ic,ic,:) = 1./Data.Indice{ic}.^2; 
                    CaractNu1(ic,ic,:) = 1./Data.Indice{ic}.^2; 
                else
                    CaractNu(ic,ic,:) = 1./(Data.Indice.^2); 
                    CaractNu1(ic,ic,:) = 1./(Data.Indice.^2); 
                end
             
        end
        MatPml = CaractEps;
end

Phys.Nu = CaractNu;
Phys.Eps = CaractEps;

for k = 1:max(Mesh.Nsd), Phys.MurMu0(1:3,1:3,k) = eye(3)*Phys.Mu0; end
for k = 1:max(Mesh.Nsd) 
    if iscell(Data.Indice)
        for ic = 1:3, Phys.EpsrEps0(ic,ic,k) = (Data.Indice{ic}(k).^2)*Phys.Eps0; end 
    else
        Phys.EpsrEps0(1:3,1:3,k) = eye(3)*(Data.Indice(k).^2)*Phys.Eps0; 
    end
end

%
CoefPml = Phys.CoefPml;
MatPmlX = [1/CoefPml 0 0 ; 0 CoefPml 0 ; 0 0 CoefPml];
MatPmlY = [CoefPml 0 0 ; 0 1/CoefPml 0 ; 0 0 CoefPml];

for k = 1:length(Phys.PmlX) 
    CaractNu(:,:,Phys.PmlX(k)) = CaractNu(:,:,Phys.PmlX(k)) * diag(1./diag(MatPmlX));
    CaractEps(:,:,Phys.PmlX(k)) = CaractEps(:,:,Phys.PmlX(k)) * MatPmlX;
    MatPml(2,2,Phys.PmlX(k)) = CoefPml;
    Phys.EpsrEps0(:,:,Phys.PmlX(k)) = Phys.EpsrEps0(:,:,Phys.PmlX(k))*MatPmlX;
    Phys.MurMu0(:,:,Phys.PmlX(k)) = Phys.MurMu0(:,:,Phys.PmlX(k))*MatPmlX;
end

for k = 1:length(Phys.PmlY) 
    CaractNu(:,:,Phys.PmlY(k)) = CaractNu(:,:,Phys.PmlY(k)) * diag(1./diag(MatPmlY));
    CaractEps(:,:,Phys.PmlY(k)) = CaractEps(:,:,Phys.PmlY(k)) * MatPmlY;
    MatPml(1,1,Phys.PmlY(k)) = CoefPml;
    Phys.EpsrEps0(:,:,Phys.PmlY(k)) = Phys.EpsrEps0(:,:,Phys.PmlY(k))*MatPmlY;
    Phys.MurMu0(:,:,Phys.PmlY(k)) = Phys.MurMu0(:,:,Phys.PmlY(k))*MatPmlY;
end

CaractNu1(1,1,:) = CaractNu(2,2,:);
CaractNu1(2,2,:) = CaractNu(1,1,:);

Phys.CaractNu1 = CaractNu1;
Phys.CaractNu = CaractNu;
Phys.CaractEps = CaractEps;
Phys.MatPml = MatPml;
if size(Mesh.CoorN,2)>1,
    [Phys.fx,Phys.fy,Phys.fz] = CoefPml2D(Mesh,Data);
end

%

if Data.ChampInc > 0
    Phys.Kx = +Phys.K0*Data.nh*sin(Data.Theta0);
    Phys.Ky = -Phys.K0*Data.nh*cos(Data.Theta0);
else
    Phys.Kx = +Phys.K0*Data.nb*sin(Data.Theta0);
    Phys.Ky = +Phys.K0*Data.nb*cos(Data.Theta0);
end


end

%%
function Phys = CaractMat3D(Mesh,Data)

% CaractMat3D 
%   Description des matériaux et prise en compte des PML
%   pour les problèmes 3D et 2D+1D solveur modal
% Formulation en E
%
% Syntaxe
%   Phys = CaractMat3D(Mesh,Data);
%
% Description
%   Mesh : structure données du maillage (CoorN, Cn, CoorA, Ca, ...,Nsd )
%   Data : données du problème
%
%   Phys : structure caractéristiques des matériaux
%       Indice : indice de réfraction des matériaux
%       PmlX : Numéros des sous-domaines de PML en x
%       PmlY : Numéros des sous-domaines de PML en y
%       PmlZ : Numéros des sous-domaines de PML en z

%   Date de la dernière version : 07 Décembre 2009 ... 23 mars 2010
%                               19 janvier 2011 (Pml non symétrique)
%                                6/2/2012 (TypePol)
%   Auteur : Mondher Besbes (LCF / CNRS / IOGS)

if size(Mesh.CoorN,2) == 3
    if numel(Data.a) == 1
        Phys.PmlX = unique(Mesh.Nsd(abs(Mesh.CoorV(:,1)) > Data.a));
    else
        Phys.PmlX = unique(Mesh.Nsd(Mesh.CoorV(:,1) < Data.a(1) | Mesh.CoorV(:,1) > Data.a(2)));
    end
    %
    if numel(Data.b) == 1
        Phys.PmlY = unique(Mesh.Nsd(abs(Mesh.CoorV(:,2)) > Data.b));
    else
        Phys.PmlY = unique(Mesh.Nsd(Mesh.CoorV(:,2) < Data.b(1) | Mesh.CoorV(:,2) > Data.b(2)));
    end
    %
    if numel(Data.c) == 1
        Phys.PmlZ = unique(Mesh.Nsd(abs(Mesh.CoorV(:,3)) > Data.c));
    else
        Phys.PmlZ = unique(Mesh.Nsd(Mesh.CoorV(:,3) < Data.c(1) | Mesh.CoorV(:,3) > Data.c(2)));
    end
elseif size(Mesh.CoorN,2) == 2
    if numel(Data.a) == 1
        Phys.PmlX = unique(Mesh.Nsd(abs(Mesh.CoorF(:,1)) > Data.a));
    else
        Phys.PmlX = unique(Mesh.Nsd(Mesh.CoorF(:,1) < Data.a(1) | Mesh.CoorF(:,1) > Data.a(2)));
    end
%
    if numel(Data.b) == 1
        Phys.PmlY = unique(Mesh.Nsd(abs(Mesh.CoorF(:,2)) > Data.b));
    else
        Phys.PmlY = unique(Mesh.Nsd(Mesh.CoorF(:,2) < Data.b(1) | Mesh.CoorF(:,2) > Data.b(2)));
    end
    Phys.PmlZ = [];
else
    disp('dimension du problème différente de 2 et de 3')
end

%
[Indice,TabPmlX,TabPmlY,TabPmlZ] = deal(Data.Indice,Phys.PmlX,Phys.PmlY,Phys.PmlZ);
Phys.TypePol = Data.TypePol;
Phys.CoefPml = Data.CoefPml;
Phys.Lambda0 = Data.Lambda0;
Phys.Indice = Data.Indice;
Phys.Sym = Data.Sym;
if isfield(Data,'SigmaG'), Phys.SigmaG = Data.SigmaG; end
if isfield(Data,'Sigma'), Phys.Sigma = Data.Sigma; end
%
Phys.c0 = 2.99792458e+14;        % µm/s
Phys.Mu0 = 4*pi*1e-7;            % Perméabilité du vide
Phys.Eps0 = 1/(36*pi*1e9);       % Permittivité du vide
Phys.Z0 = sqrt(Phys.Mu0/Phys.Eps0);
Phys.K0 = 2*pi/Phys.Lambda0;

CaractNu = zeros(3,3,max(Mesh.Nsd));
CaractEps = zeros(3,3,max(Mesh.Nsd));

for k = 1:max(Mesh.Nsd), Phys.Mu(1:3,1:3,k) = eye(3)*Phys.Mu0; end
%
for k = 1:max(Mesh.Nsd), CaractNu(1:3,1:3,k) = eye(3); end
%
for ic = 1:3
    if iscell(Data.Indice)
        CaractEps(ic,ic,:) = Indice{ic}.^2; 
        Phys.Eps(ic,ic,:) = Indice{ic}.^2*Phys.Eps0;
    else
        CaractEps(ic,ic,:) = Indice.^2;
        Phys.Eps(ic,ic,:) = Indice.^2*Phys.Eps0;
    end
end

for k = 1:max(Mesh.Nsd), Phys.MurMu0(1:3,1:3,k) = eye(3)*Phys.Mu0; end
for k = 1:max(Mesh.Nsd) 
    if iscell(Data.Indice)
        for ic = 1:3, Phys.EpsrEps0(ic,ic,k) = (Data.Indice{ic}(k).^2)*Phys.Eps0; end 
    else
        Phys.EpsrEps0(1:3,1:3,k) = eye(3)*(Data.Indice(k).^2)*Phys.Eps0; 
    end
end

CoefPml = Phys.CoefPml;
MatPmlX = [1/CoefPml 0 0 ; 0 CoefPml 0 ; 0 0 CoefPml];
MatPmlY = [CoefPml 0 0 ; 0 1/CoefPml 0 ; 0 0 CoefPml];
MatPmlZ = [CoefPml 0 0 ; 0 CoefPml 0 ; 0 0 1/CoefPml];

for k = 1:length(TabPmlX) 
    CaractNu(:,:,TabPmlX(k)) = CaractNu(:,:,TabPmlX(k)) * diag(1./diag(MatPmlX));
    CaractEps(:,:,TabPmlX(k)) = CaractEps(:,:,TabPmlX(k)) * MatPmlX;
    Phys.EpsrEps0(:,:,TabPmlX(k)) = Phys.EpsrEps0(:,:,TabPmlX(k))*MatPmlX;
    Phys.MurMu0(:,:,TabPmlX(k)) = Phys.MurMu0(:,:,TabPmlX(k))*MatPmlX;
end

for k = 1:length(TabPmlY) 
    CaractNu(:,:,TabPmlY(k)) = CaractNu(:,:,TabPmlY(k)) * diag(1./diag(MatPmlY));
    CaractEps(:,:,TabPmlY(k)) = CaractEps(:,:,TabPmlY(k)) * MatPmlY;
        Phys.EpsrEps0(:,:,TabPmlY(k)) = Phys.EpsrEps0(:,:,TabPmlY(k))*MatPmlY;
    Phys.MurMu0(:,:,TabPmlY(k)) = Phys.MurMu0(:,:,TabPmlY(k))*MatPmlY;

end

for k = 1:length(TabPmlZ) 
    CaractNu(:,:,TabPmlZ(k)) = CaractNu(:,:,TabPmlZ(k)) * diag(1./diag(MatPmlZ));
    CaractEps(:,:,TabPmlZ(k)) = CaractEps(:,:,TabPmlZ(k)) * MatPmlZ;
    Phys.EpsrEps0(:,:,TabPmlZ(k)) = Phys.EpsrEps0(:,:,TabPmlZ(k))*MatPmlZ;
    Phys.MurMu0(:,:,TabPmlZ(k)) = Phys.MurMu0(:,:,TabPmlZ(k))*MatPmlZ;

end

Phys.CaractNu = CaractNu;
Phys.CaractEps = CaractEps;

% Cas du solveur modal 2D
CaractNu1 = zeros(size(CaractNu));
if size(Mesh.CoorN,2) == 2
    for k = 1:max(Mesh.Nsd) 
        [CaractNu1(1,1,k),CaractNu1(2,2,k)] = ...
            deal(CaractNu(2,2,k),CaractNu(1,1,k)); 
    end
    Phys.CaractNu1 = CaractNu1;
end


if size(Mesh.CoorN,2) == 3
    if Data.ChampInc > 0
        Phys.Kx = +Phys.K0*Data.nh*sin(Data.Theta0)*cos(Data.Phi0);
        Phys.Ky = +Phys.K0*Data.nh*sin(Data.Theta0)*sin(Data.Phi0);
        Phys.Kz = -Phys.K0*Data.nh*cos(Data.Theta0);
    else
        Phys.Kx = +Phys.K0*Data.nb*sin(Data.Theta0)*cos(Data.Phi0);
        Phys.Ky = +Phys.K0*Data.nb*sin(Data.Theta0)*sin(Data.Phi0);
        Phys.Kz = +Phys.K0*Data.nb*cos(Data.Theta0);
    end
elseif size(Mesh.CoorN,2) == 2
    if Data.ChampInc > 0
        Phys.Kx = +Phys.K0*Data.nh*sin(Data.Theta0);
        Phys.Ky = -Phys.K0*Data.nh*cos(Data.Theta0);
    else
        Phys.Kx = +Phys.K0*Data.nb*sin(Data.Theta0);
        Phys.Ky = +Phys.K0*Data.nb*cos(Data.Theta0);
    end
end

end
