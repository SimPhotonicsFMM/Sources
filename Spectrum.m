function varargout = Spectrum(index,geom,lambda,theta,inc,varargin)

% Spectrum 
%   Parallel calculation of the spectrum of a photonic structure  for 
%   the 2 polarizations TM and TE illuminated from above or below 
%   - Fourier Modal Method -
%
% Syntax
%   Model = Spectrum(index,geom,lambda,theta,inc)
%   [R_tem,T_tem] = Spectrum(index,geom,lambda,theta,inc)
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc)
%   varargout = Spectrum(index,geom,lambda,theta,inc,'param1',val1,...)
%
% Description
%   index : Refractive indices from top to bottom [nh nc_1...nc_i nb]
%   geom :  Thickness of different layers from top to bottom [hc_1..hc_i..]
%        or structure geometric data (dx,dy,hc,mn,ab,Angle,Dep,Np)
%        with superformula or for rectangular inclusions (dx,dy,hc,lix,liy)
%
%                     nh
%               ---------------
%                   nc_1,hc_1
%               ---------------
%                    ...
%                  nc_i,hc_i
%                    ...
%               ---------------
%                     nb
%
%   lambda : Wavelength (µm)
%   theta  : Incidence angle (rd)
%   inc    : Incidence from top =+1 , from bottom =-1
%
%   Additional parameters: 'param1',val1,... 
%       mx  : Number of Fourier terms in x, by default mx=0
%       my  : Number of Fourier terms in y, by default my=0
%       Nper: Possible number of periods (Bragg miror), by default Nper=1
%       Num : Diffracted order numbers, by default Num=[]
%       Phi0 : Azimuthal angle (rd), by default Phi0=0
%       SymY : Symmetry in y, =0 :TM (PMC), =1 : TE (PEC), by default =2
%              (SymY=0 or 1: Only total R/T will be computed)
%       Nsub : Number of subdivisions for FD_FMM calculation
%
%       dx  : Period in x, , by default dx=lambda(1)
%       dy  : Period in y, , by default dy=lambda(1)
%       lix : Width of inclusions in x, , by default lix=[]
%       liy : Width of inclusions in y, , by default liy=[]
%
%       nx  : Number of nodes in x, , by default nx=2
%       ny  : Number of nodes in y, , by default ny=2
%       nz  : Number of nodes in z, , by default nz=2
%
%   R_tem: Array of total reflectivity TM-TM, TM-TE, TE-TM, TE-TE
%   T_tem: Array of total transmittivity TM-TM, TM-TE, TE-TM, TE-TE
%   R_tm : Array of total reflectivity TM (or diffraction order with 'Num') 
%   R_te : Array of total reflectivity TE
%   T_tm : Array of total transmittivity TM
%   T_te : Array of total transmittivity TE
%   Model: = R_tem if length(Theta or/and lambda)>1 
%          otherwize a Matlab structure with fields 
%          (Mesh,Phys,Data,R_tem,T_tem,Sb,Sh,MatS,CoefD)
%
% Examples : 
%   Dioptre
%   -------
%   index = [1 1.5]; geom = []; lambda = 0.6; theta = (0:1:89.9)*pi/180; inc = -1; 
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc);
%   figure, PlotCoefRTA(lambda,theta,R_tm,R_te), legend('R_{tm} (%)','R_{te} (%)')
%
%   AntiReflection Coating
%   ----------------------
%   index = [1 2 4]; geom = 0.08; lambda = .4:.01:.8; theta = pi/12; inc = 1; 
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc);
%   figure, PlotCoefRTA(lambda,theta,R_tm,T_tm), legend('R_{tm} (%)','T_{tm} (%)')
%
%   Bragg Miror
%   -----------
%   ld0 = .6; nc1 = 1.5; nc2 = 2.2; hc1 = ld0/4/nc1; hc2 = ld0/4/nc2; Np = 20;
%   index = [1 nc1 nc2 1.5]; geom = [hc1 hc2]; lambda = .4:.001:.8; theta = 0; inc = 1; 
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc,'Nper',Np);
%   figure, PlotCoefRTA(lambda,theta,R_tm,T_tm), legend('R_{tm} (%)','T_{tm} (%)')
%
%   SPR Biosensor
%   -------------
%   index = {1.33 , @(x) IndexVal('Au',x) , IndexVal('SF10')}; 
%   geom = 0.05; lambda = .85; theta = linspace(50,60,201)*pi/180; inc = -1; 
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc);
%   figure, PlotCoefRTA(lambda,theta,R_tm,T_tm), legend('R_{tm} (%)','T_{tm} (%)')
%
%   Grating SPR biosensor 1D ou 2D
%   ------------------------------
%   lambda = 0.85; %linspace(.7,.9,41); 
%   theta = linspace(50,60,51)*pi/180;
%   index = {1.33 , {@(x) IndexVal('Au',x) 1.33} ,@(x) IndexVal('Au',x), 1.7};
%   geom = [.015 0.04];  inc = -1; 
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc,...
%                        'dx',.2,'dy',.2,'lix',{0.1 0},'liy',{0 0},'mx',10,'my',0);
%   %geom = SetGeom('dx',.2,'dy',.2,'hc',[.015 .04],'mn',{4, 4},'ab',{[.05 .1] [.1 .1]});
%   %geom = SetGeom('dx',.2,'dy',.2,'hc',[.015 .04],'lix',{.1 , [] },'liy',{.2 []});
%   %[R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc,'mx',10,'my',0);
%   figure, PlotCoefRTA(lambda,theta,R_tm,T_tm), legend('R_{tm} (%)','T_{tm} (%)')
%
%   MIM plasmonic antenna
%   ---------------------
%   lambda = linspace(2.5,4.5,61); theta = 0; inc = +1; 
%   index = {1 , {@(x) IndexVal('Au',x) 1} , IndexVal('SiO2'), @(x) IndexVal('Au',x)};
%   geom = SetGeom('dx',2,'dy',2,'hc',[.05 0.22],'lix',{1, []},'liy',{.1 []});
%   [R_tm,T_tm,R_te,T_te] = Spectrum(index,geom,lambda,theta,inc,'mx',7,'my',7);
%   figure, PlotCoefRTA(lambda,theta,R_tm,1-R_tm-T_tm), legend('R_{tm} (%)','A_{tm} (%)')

% Date of the latest version : 17 Avril 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

%
% Initialiser la structure des paramètres
global Index SaveOption

if iscell(index)
    for k=1:length(index) 
        if ~isnumeric(index{k}) 
            load('IndexData.mat','Index')   % Index structure
            break
        end
    end
end
lambda = lambda + 0*eps;

%
s = struct('Nper',1,'mx',0,'my',0,'Phi0',0,'npx',2,'npy',2,'npz',2,'Num',[],...
     'dx',lambda(1)/3,'dy',lambda(1)/3,'lix',[],'liy',[],'Sym',[2 2],'SymY',2,...
     'Nsub',[],'Eps',1e-6);
%
% Listes des paramètres
%if nargin == 6, s.Nper = varargin{1}; end
if numel(varargin) == 1, varin = varargin{1}; else, varin = varargin; end
if nargin >= 6
    f = varin(1:2:end);
    v = varin(2:2:end);
    FieldData = fieldnames(s);
    for k = 1:length(f)
        P = ismember(FieldData,f(k));
        if any(P) 
            s.(FieldData{P}) = v{k};
        else % nouveau champ
            s.(f{k}) = v{k};
            %warning(['Nouvelle donnée : ', f{k}]); 
        end
    end
end
Nper = s.Nper;
%
% Indices des milieux haut et bas et des differentes couches
if iscell(index)
    nh = index{1};     % Indice milieu haut
    nb = index{end};   % Indice milieu bas
    nc = index(end-1:-1:2);
else
    nh = index(1);     % Indice milieu haut
    nb = index(end);   % Indice milieu bas
    nc = index(end-1:-1:2);
end
%
% Hauteurs des couches et largeurs des inclusions
[dx,dy] = deal(Inf);
if isstruct(geom) && isfield(geom,'Cn')
    Mesh = geom;
else
if isstruct(geom)
%     if isfield(geom,'npx')
%         [s.npx,s.npy,s.npz] = deal(geom.npx,geom.npy,geom.npz);
%     end
    %
    if isfield(geom,'Xpml')
        dx = geom.dx;
        geom.dx = dx+geom.Xpml;
        dy = geom.dy;
        geom.dy = dy+geom.Ypml;
    end
    %
    Mesh = MeshLayer(geom);
    %
    if isfield(geom,'Xpml')
        Mesh.Nsd(abs(Mesh.CoorV(:,2))>dy/2) = max(Mesh.Nsd)+1;
        Mesh.Nsd(abs(Mesh.CoorV(:,2))>dy/2 & abs(Mesh.CoorV(:,1))<dx/2) = max(Mesh.Nsd)+1;
    end
else
    hc = geom(end:-1:1);
    if isempty(nc), nc = nb; hc = eps; end    % Cas d'un dioptre
    %
    Nb_hc = length(hc);
    if Nb_hc>1 && isempty(s.lix), s.lix = cell(Nb_hc,1); end
    if Nb_hc>1 && isempty(s.liy), s.liy = cell(Nb_hc,1); end
    % Description de la géométrie
    SaveOption = 1;
    geom = []; geom.dx = s.dx; geom.dy = s.dy;
    Mesh = MeshLayer(s.dx,s.lix(end:-1:1),s.dy,s.liy(end:-1:1),hc(:)',...
                     s.npx(end:-1:1),s.npy(end:-1:1),s.npz(end:-1:1));
end
end
%
% Description des données
if iscell(nc), Indice = nc; else, Indice = num2cell(nc); end
if isfield(s,'SymY'), s.Sym = [2 s.SymY]; end
%
Data = SetData('Lambda0',lambda(1),'Theta0',theta(1),'Phi0',s.Phi0(1),...
       'ChampInc',inc,'mx',s.mx,'my',s.my,'nh',nh,'nb',nb,'Indice',Indice,...
       'a',[-dx/2 dx/2],'b',[-dy/2 dy/2],'CoefPml',1+1i,'Sym',s.Sym,'Eps',s.Eps);
%
if Data(1).Sym(2) ~= 2 && Data(1).Phi0 ~= 0
    error('if Phi0 is not equal 0, there is no symmetry in y')
end

%
if Data(1).Sym(2) ~= 2 && ~isempty(s.Num)
    error('if Num is not empty, symmetry is not authorised ')
end

%
if ~isempty(s.Num), NbOrder = length(s.Num); else, NbOrder = 1; end
if isfield(s,'Field') && ~isempty(s.Field), for k=1:length(Data), Data(k).Field = s.Field; end, end
 

NumLayer = [];
if isfield(s,'Nsub') && ~isempty(s.Nsub)
    if numel(find(s.Nsub~=0))>1 || (numel(find(s.Nsub~=0))==1 && Nper>1) 
        %warning('FD_FMM is only used for one layer'); 
    end
    for k=1:length(Data) 
        Data(k).Nsub = s.Nsub(end-k+1); 
        if Data(k).Nsub ~= 0, 
            NumLayer = [NumLayer k];
            if isfield(s,'POD'), Data(k).POD = s.POD; end
        end
    end 

end
%


% Calcul du spectre en fonction de Lambda et Theta
[Xl,Xa] = ndgrid(lambda,theta);

NbArgOut = nargout;
if NbArgOut == 4
    [R_tm,T_tm,R_te,T_te] = deal(zeros(length(Xl(:)),NbOrder));
else
    if sum(Data(1).Sym) ~= 4
        [R_tem,T_tem] = deal(zeros(length(Xl(:)),1));
    else
        [R_tem,T_tem] = deal(zeros(length(Xl(:)),4*NbOrder));
    end
end
%
% Check if h == lambda/4n
for kh = 1:length(Mesh)
    hc = max(Mesh(kh).CoorN(:,end))-min(Mesh(kh).CoorN(:,end));
    for k = 1:length(lambda)  % Loop lambda
        Data1 = InterpIndex(Data(kh),lambda(k)); % Mise à jour des indices, nb et nh
        if isscalar(Data1.Indice) && any(abs(theta) < pi/36)
            h0 = lambda(k)/(4*Data1.Indice);
            if abs(hc/h0 - round(hc/h0))<=eps
                if mod(hc/h0,2) ~= 0
                    error(['Divide into 2 slices the film of thickness = ' num2str(hc)])
                end
            end
        end
    end
end


%
if isscalar(lambda) && isscalar(theta)
    Data = SetData(Data,'Lambda0',lambda,'Theta0',theta);
    Data = InterpIndex(Data); % Mise à jour des indices, nb et nh

    %
    Phys = CaractMat(Mesh,Data);
    if isfield(s,'Kx') && ~isfield(s,'Ky'), Phys = CaractMat(Mesh,Data,s.Kx,0); end
    if ~isfield(s,'Kx') && isfield(s,'Ky'), Phys = CaractMat(Mesh,Data,0,s.Ky); end
    if isfield(s,'Kx') && isfield(s,'Ky'), Phys = CaractMat(Mesh,Data,s.Kx,s.Ky); end
    Sb = CalculMatS(Data,Mesh,Phys,-1); % milieu bas
    Sh = CalculMatS(Data,Mesh,Phys,+1); % milieu haut
    MatS = CalculMatS(Data,Mesh,Phys); % Matrices S des différentes couches
    if Nper > 1, MatS0 = ProdMatS(MatS,Nper); else, MatS0 = MatS; end
    if isempty(Phys(1).PmlX) && isempty(Phys(1).PmlY)
        if ~isempty(NumLayer)
            [r,t,CoefD,R,T,E,H] = CalculCoefRT(Sb,MatS0,Sh,NumLayer);
            %n = length(MatS{1,1});
            if length(NumLayer) == 1
                MatS{NumLayer,6} = [E H];
            else
                P = 0;
                for k = 1:length(NumLayer)
                    Nsub = Data(NumLayer(k)).Nsub;
                    E0 = []; H0 = [];
                    for ks = 1:Nsub+1
                        E0 = [E0 ; E{ks+P}]; 
                        H0 = [H0 ; H{ks+P}]; 
                    end
                    MatS{NumLayer(k),6} = [E0 H0];
                    P = P + Nsub;
                end
            end
        else
            [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS0,Sh);
        end
        if sum(Data(1).Sym) ~= 4
            [R_tem,T_tem] = CoefRTA(r,t);
        else
            if ~isempty(s.Num)
                [R_tem,T_tem] = deal(zeros(length(s.Num),4));
                for kn = 1:length(s.Num)
                    [R_tem(kn,:),T_tem(kn,:)] = CoefRTA(r,t,s.Num(kn));
                end
            else
                [R_tem,T_tem] = CoefRTA(r,t);
            end
        end
    else
        CoefD = [];
    end

else
    parfor k = 1:length(Xl(:))  % Loop lambda & theta
        %
        Data0 = SetData(Data,'Lambda0',Xl(k),'Theta0',Xa(k),'nb',nb,'nh',nh);
        Data1 = InterpIndex(Data0); % Mise à jour des indices, nb et nh
        %
        %
        Phys1 = CaractMat(Mesh,Data1);
        [r,t] = CalculFMM(Data1,Mesh,Phys1,Nper);
        if sum(Data1(1).Sym) ~= 4
            [R,T] = CoefRTA(r,t);
            if NbArgOut == 4
                if Data(1).Sym(2) == 0, [R_tm(k,:),T_tm(k,:)] = deal(R(:,1),T(:,1)); end
                if Data(1).Sym(2) == 1, [R_te(k,:),T_te(k,:)] = deal(R(:,1),T(:,1)); end
            elseif NbArgOut <= 2
                 R_tem(k,:) = R(:,1);
                 T_tem(k,:) = T(:,1);
            end
        else
            if ~isempty(s.Num)
                [R,T] = deal(zeros(length(s.Num),4));
                for kn = 1:length(s.Num)
                    [R(kn,:),T(kn,:)] = CoefRTA(r,t,s.Num(kn));
                end
            else
                [R,T] = CoefRTA(r,t);
            end
            if NbArgOut == 4
                [R_tm(k,:),T_tm(k,:)] = deal(R(:,1)+R(:,2),T(:,1)+T(:,2));
                [R_te(k,:),T_te(k,:)] = deal(R(:,3)+R(:,4),T(:,3)+T(:,4));
            elseif NbArgOut <= 2
                R_tem(k,:) = R(:);
                T_tem(k,:) = T(:);
            end
        end
    end
end

if NbArgOut == 1
    if numel(lambda) == 1 && numel(theta) == 1
        Model = struct('Data',Data,'Mesh',Mesh,'Phys',Phys,...
                        'R_tem',R_tem,'T_tem',T_tem,'CoefD',CoefD,...
                        'R',R,'T',T);
        Model.Sb = Sb; Model.Sh = Sh; Model.MatS = MatS;
        if Nper>1, Model.Nper = Nper; end
    else
        Model = R_tem;
    end
    varargout = {Model};
elseif NbArgOut == 4
        
    if numel(lambda) == 1 && numel(theta) == 1
        if sum(Data(1).Sym) ~= 4
            [R_tm,T_tm] = deal(R_tem,T_tem);
            [R_te,T_te] = deal(R_tem,T_tem);
        else
            [R_tm,T_tm] = deal(R_tem(:,1)+R_tem(:,2),T_tem(:,1)+T_tem(:,2));
            [R_te,T_te] = deal(R_tem(:,3)+R_tem(:,4),T_tem(:,3)+T_tem(:,4));
        end
    end
    %
    varargout = {squeeze(reshape(R_tm,size(Xl,1),size(Xl,2),NbOrder)), ...
                 squeeze(reshape(T_tm,size(Xl,1),size(Xl,2),NbOrder)), ...
                 squeeze(reshape(R_te,size(Xl,1),size(Xl,2),NbOrder)), ...
                 squeeze(reshape(T_te,size(Xl,1),size(Xl,2),NbOrder))};
    %varargout = {R_tm T_tm R_te T_te};
    
else
    if sum(Data(1).Sym) ~= 4
        varargout = {R_tem T_tem};
    else
    varargout = {squeeze(reshape(R_tem,size(Xl,1),size(Xl,2),NbOrder,4)), ...
                 squeeze(reshape(T_tem,size(Xl,1),size(Xl,2),NbOrder,4))};
    end
end

end

