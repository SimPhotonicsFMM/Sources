function [MatS,Sdm] = CalculMatS3D(Data,Mesh,Phys,Dof,Mat,SdMembre,Sdm0)

% CalculMatS3D 
%   Calculation of S-Matrix MMF or HYB model for 3D Geometry (2D Grating)
%       Symmetry at y-axis : Data.Sym = [2 0] TM , = [2 1] TE
%
% Syntax
%   MatS = CalculMatS3D(Data,Mesh,Phys);
%   MatS = CalculMatS3D(Data,Mesh,Phys,Coef);
%   MatS = CalculMatS3D(Data,Mesh,Phys,Dof,Mat);
%
% Description
%   Data : Data of the priblem
%   Mesh : Mesh of the structure (CoorN, Cn, Nsd,...CoorA, Ca, ExtAr, )
%   Phys : Material properties (Epsr, Mur, ...,Lambda0, K0, Omega, TypePol)
%   Coef : +1 Upstrat S-Matrix ; -1 Substrat S-Matrix
%
%   MatS : S-matrix   [Eh;Hb] = MatS*[Eb;Hh]

% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

Coef = 0;
if nargin == 4 && ~isstruct(Dof), 
    Coef = Dof; 
    Data = Data(1);
    Mesh = Mesh(1);
    Phys = Phys(1);
end
if isfield(Data,'mx') && isfinite(Data.mx), mx = Data.mx;  else mx = 0; end
if isfield(Data,'my') && isfinite(Data.my), my = Data.my; else my = 0; end

pol = Phys.TypePol;
m = 2*mx+1;
Data.mx = mx;
Data.my = my;

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


%% Calcul des matrices S

if Coef == -1
% --------------- côté bas -----------------------------------------------
    %
    if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
        Data.Indice(:) = Data.nb;
        Phys = CaractMat(Mesh,Data);
        MatSb = CalculMatS3D(Data,Mesh,Phys);
        [P,Q,V,Vp] = deal(MatSb{3},MatSb{4},MatSb{5},MatSb{6});
        n = size(MatSb{1},1)/2;
        InvSb12 = inv(MatSb{1}(1:n,n+1:end));
    else

    Data.nc = Data.nb;
    %
    [P,Q,V,Vp] = CalculVVp(Data,Phys,BetaX,BetaY,Dx,Dy);
    InvP = inv(P);
    InvQ = inv(Q);
    %
    if sum(Data.Sym) ~= 4
        [R_E,R_H,InvR_E,InvR_H] = MatSym(Data);
        normE = sqrt(sum(R_E.^2,1));
        normH = sqrt(sum(R_H.^2,1));
        P1 = InvR_H*P*R_H*diag(1./normH);
        P2 = InvR_H*P*R_E*diag(1./normE);
        InvP1 = diag(normH)*InvR_H*InvP*R_H;
        InvP2 = diag(normE)*InvR_E*InvP*R_H;
        Q1 = InvR_E*Q*R_E*diag(1./normE);
        Q2 = InvR_E*Q*R_H*diag(1./normH);
        InvQ2 = diag(normH)*InvR_H*InvQ*R_E;
        %
%         if Data.Sym(1) ~= 2 && Data.Sym(2) == 2
%         Vp = diag(InvR_E*diag(Vp)*R_E);
%         else
        Vp = diag(InvR_H*diag(Vp)*R_H);
%        end

        %
%         MatSb = {[Q1+Q2*InvP1*P2 , -Q2*InvP1 ; InvP1*P2 , -InvP1], ...
%                  zeros(length(Q1)*2,1)};

if Data.Sym(1) == 2 && Data.Sym(2) ~= 2
        MatSb = {[2*Q2 , -Q2*InvP1 ; eye(size(InvP1*P2)) , -InvP1], ...
                 zeros(length(Q1)*2,1)};
         else
% 
         MatSb = {[2*Q1 , -Q2*InvP1 ; eye(size(InvP1*P2)) , -InvP1], ...
                  zeros(length(Q1)*2,1)};
         end
        InvSb12 = -P1*InvQ2;
    else
        MatSb = {[2*Q , -Q*InvP ; eye(size(Q)) , -InvP], ...
            zeros(length(P)*2,1)};
        InvSb12 = -P*InvQ;
    end
    %


    %
    end
    %
    Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
    Pd = Pvp(:).';
    %
    n = size(MatSb{1},1)/2;

    if Data.ChampInc == -1

        if sum(Data.Sym) ~= 4
            Pdi = find(abs(abs(Vp)-abs(Data.nb*Phys.K0*cos(Data.Theta0)))<Phys.K0/1e6)';
        else
            Pdi = [(length(mx)-1)/2+1+length(mx) (length(mx)-1)/2+1];
        end
         if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
             Pdi = find(abs(abs(Vp)-abs(Data.nb*Phys.K0))<Phys.K0/1e6)';
             Pdi = Pdi(2:-1:1);
             Pd = 1:length(Vp);
             %Pdi = Pd;
         end
%        Pdi = [(length(mx)-1)/2+1 (length(mx)-1)/2+1+length(mx)];
    else
        Pdi = [];
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY), Pd = 1:length(Vp); end
    end
    %
    Sb{1} = MatSb{1}([1:n,Pd+n],[Pdi,n+1:end]);
    Sb{2} = MatSb{2}([1:n,Pd+n]);
    [Sb{3},Sb{4},Sb{5},Sb{6},Sb{7},Sb{8},Sb{9}] = deal(P,Q,V,Vp,Pdi,Pd,InvSb12);
    MatS = Sb;

elseif Coef == +1
% --------------- côté haut ------------------------------------------------
    %
    if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
        Data.Indice(:) = Data.nb;
        Phys = CaractMat(Mesh,Data);
        MatSh = CalculMatS3D(Data,Mesh,Phys);
        [P,Q,V,Vp] = deal(MatSh{3},MatSh{4},MatSh{5},MatSh{6});
    else

    Data.nc = Data.nh;
    Data.Alfa = [];

    [P,Q,V,Vp] = CalculVVp(Data,Phys,BetaX,BetaY,Dx,Dy);

    InvQ = inv(Q);
    if sum(Data.Sym) ~= 4
        InvP = inv(P);
        [R_E,R_H,InvR_E,InvR_H] = MatSym(Data);
        normE = sqrt(sum(R_E.^2,1));
        normH = sqrt(sum(R_H.^2,1));
        P1 = InvR_H*P*R_H*diag(1./normH);
        P2 = InvR_H*P*R_E*diag(1./normE);
        %Q2 = InvR_E*Q*R_H*diag(1./normH);
        InvQ1 = diag(normE)*InvR_E*InvQ*R_E;
        InvQ2 = diag(normH)*InvR_H*InvQ*R_E;
        InvP2 = diag(normE)*InvR_E*InvP*R_H;
%         if Data.ChampInc == -1
%             Vp = diag(InvR_E*diag(Vp)*R_E);
%         elseif Data.ChampInc == +1
            Vp = diag(InvR_H*diag(Vp)*R_H);
%        end
        %
%         MatSh = {[InvQ1 , -InvQ1*Q2 ; P2*InvQ1 , -(P2*InvQ1*Q2+P1)], ...
%                  zeros(length(P1)*2,1)};
%         MatSh = {[InvQ1/2+InvP2*P1*InvQ2/2 , -InvP2*P1 ; P1*InvQ2 , -2*P1], ...
%                  zeros(length(P1)*2,1)};

if Data.Sym(1) == 2 && Data.Sym(2) ~= 2           
            MatSh = {[InvQ2 , -eye(size(InvP2*P1)) ; P1*InvQ2 , -2*P1], ...
                 zeros(length(P1)*2,1)};
         else
             MatSh = {[InvQ1 , -eye(size(InvP2*P1)) ; P1*InvQ2 , -2*P1], ...
                  zeros(length(P1)*2,1)};
% 
         end
    else
        MatSh = {[InvQ , -eye(size(Q)); P*InvQ , -2*P], ...
                 zeros(length(P)*2,1)};
    end

    %


    end
    %% 
    Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
    Pu = Pvp.';

    n = size(MatSh{1},1)/2;
    if Data.ChampInc == +1
        if sum(Data.Sym) ~= 4
            Pui = find(abs(abs(Vp)-abs(Data.nh*Phys.K0*cos(Data.Theta0)))<Phys.K0/1e6)'+n;

        else
            Pui = [(length(mx)-1)/2+1+length(mx)+n (length(mx)-1)/2+1+n];
            if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
                Pui = find(abs(abs(Vp)-abs(Data.nh*Phys.K0))<Phys.K0/1e6)'+n;
                Pui = Pui(2:-1:1);
                Pu = 1:length(Vp);
                %Pui = Pu;
            end
        end
%        Pui = [(length(mx)-1)/2+1+n (length(mx)-1)/2+1+length(mx)+n]; 
    else
        Pui = [];
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY), Pu = 1:length(Vp); end
    end
    %
    Sh{1} = MatSh{1}([Pu,n+1:2*n],[1:n,Pui]);
    Sh{2} = MatSh{2}([Pu,n+1:2*n]);
    [Sh{3},Sh{4},Sh{5},Sh{6},Sh{7},Sh{8}] = deal(P,Q,V,Vp,Pui,Pu);
    MatS = Sh;

else
% --------------- couche homogène ou réseau ------------------------------------------------


    if nargin == 3     % MMF
        %
        [Mesh.Cn,Mesh.CoorN,Mesh.Nsd,Mesh.CoorV]=deal(Mesh.Cn,Mesh.CoorN,Mesh.Nsd,Mesh.CoorV);
        Mesh = ShiftMesh(utilMesh2(Mesh));
        %
        if numel(Data.Indice) == 1 
            Data.nc = Data.Indice;
            if sum(Data.Sym) ~= 4, Data.lx = 1; Data.ly = 1; end
        else
            Pnx = abs(Mesh.CoorN(:,2)-min(Mesh.CoorN(:,2)))<Dx/1e6 & abs(Mesh.CoorN(:,3)-min(Mesh.CoorN(:,3)))<Dx/1e6;
            x = sort(Mesh.CoorN(Pnx,1)-min(Mesh.CoorN(:,1)));
            Data.lx = unique(x(2:end))/Dx;
            Pny = abs(Mesh.CoorN(:,1)-min(Mesh.CoorN(:,1)))<Dx/1e6 & abs(Mesh.CoorN(:,3)-min(Mesh.CoorN(:,3)))<Dx/1e6;
            y = sort(Mesh.CoorN(Pny,2)-min(Mesh.CoorN(:,2)));
            Data.ly = unique(y(2:end))/Dy;
            %
            Pn = find(abs(Mesh.CoorN(:,3)-min(Mesh.CoorN(:,3)))<Dx/1e6);
            Pe = unique(find(sum(ismember(Mesh.Cn,Pn),2)));
            [~,P]=sort(Mesh.CoorV(Pe,2));
            nc = Data.Indice(Mesh.Nsd(Pe(P))); 
            Data.nc = reshape(nc(:),length(Data.lx),length(Data.ly));
            %
            if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
                for k = 1:3
                    [Pk,Qk] = deal(squeeze(Phys.CaractEps(k,k,:)),1./squeeze(Phys.CaractNu(k,k,:)));
                    [Pk,Qk] = deal(Pk(Mesh.Nsd(Pe(P))),Qk(Mesh.Nsd(Pe(P))));
                    Data.MatP{k} = reshape(Pk(:),length(Data.lx),length(Data.ly));
                    Data.MatQ{k} = reshape(Qk(:),length(Data.lx),length(Data.ly));
                end
            end
        end
        %
        [P,Q,V,Vp,InvEpz,InvMuz] = CalculVVp(Data,Phys,BetaX,BetaY,Dx,Dy);
        %
        hc = max(Mesh.CoorN(:,3))-min(Mesh.CoorN(:,3)); %Data.hc;

        if isfield(Data,'Nsub') && ~isempty(Data.Nsub) && Data.Nsub ~= 0
            Data.hc = hc;
            [MatS{1},MatS{2},MatS{3},MatS{4},MatS{5},MatS{6},MatS{7}] = deal(P,Q,Data,Phys,[],[],[InvEpz(:),InvMuz(:)]);
            return
        end


%        InvP = inv(P);
%        InvQ = inv(Q);
        try
            dP = decomposition(P,'lu');
            dQ = decomposition(Q,'lu');
        catch
            dP = P;
            dQ = Q;
        end

         %
%         MatS{1} = [Q*diag(1./cosh(Vp*hc))*InvQ , Q*diag(tanh(Vp*hc))*InvP;
%                   -P*diag(tanh(Vp*hc))*InvQ    , P*diag(1./cosh(Vp*hc))*InvP];
        %
%         InvCosH = diag(1./cosh(Vp*hc));
%         TanH = diag(tanh(Vp*hc));
%         MatS{1} = [(Q*InvCosH)*InvQ , (Q*TanH)*InvP;
%                   -(P*TanH)*InvQ    , (P*InvCosH)*InvP];
        %
        InvCosH = repmat(reshape(1./cosh(Vp*hc),1,length(Vp)),length(Vp),1);
        TanH = repmat(reshape(tanh(Vp*hc),1,length(Vp)),length(Vp),1);
%         MatS{1} = [(Q.*InvCosH)*InvQ , (Q.*TanH)*InvP;
%                   -(P.*TanH)*InvQ    , (P.*InvCosH)*InvP];
        MatS{1} = [(Q.*InvCosH)/dQ , (Q.*TanH)/dP;
                  -(P.*TanH)/dQ    , (P.*InvCosH)/dP];

        MatS{2} = zeros(length(P)*2,1);
        [MatS{3},MatS{4},MatS{5},MatS{6},MatS{7}] = deal(P,Q,V,Vp,[InvEpz(:),InvMuz(:)]);
    else                % HYB
        if nargin == 4
            tic, [K,M,Pe0] = InitCalculMat(Data,Mesh,Dof); toc
            Mat = CalculMat(Mesh,Dof,Phys,K,M,Pe0);
        end

        if nargin == 7
            [Eh,Hh,Eb,Hb] = deal(Sdm0.Eh,Sdm0.Hh,Sdm0.Eb,Sdm0.Hb);
        else
            [Eh,Hh,Eb,Hb] = CoefFourier2D(Data,Mesh,Phys,Dof);
        end
        %
        if iscell(Mat)
            [K,M,Pe0] = deal(Mat{1},Mat{2},Mat{3});
            Mat = CalculMat(Mesh,Dof,Phys,K,M,Pe0);
        end
        %
        nb = size(Eb,2);
        nh = size(Eh,2);
        n = nb+nh;
        n0 = length(Mat)-n;
        %
        Zero = zeros(m);
        [L,U,M1,M2] = lu(Mat(1:n0,1:n0));
        if nargin == 7
            MatM = Sdm0.MatM;
        else
            MatM = inv(Mat(n0+1:n0+n,n0+1:n0+n) - ((Mat(n0+1:n0+n,1:n0)*M2/U)...
                       *(L\(M1*Mat(1:n0,n0+1:n0+n)))));
        end
        %           
%         MatG.bas = [Zero,   -Eh*MatM(1:nh,nh+1:nh+nb)*Hb;...
%                     -eye(m),  -Eb*MatM(nh+1:nh+nb,nh+1:nh+nb)*Hb];
%         %        
%         MatG.haut = [eye(m), -Eh*MatM(1:nh,1:nh)*Hh; ...
%                      Zero,   -Eb*MatM(nh+1:nh+nb,1:nh)*Hh];
        %
        if nargin < 6,
            if isfield(Data,'DofType') && strcmp(Data.DofType,'Edge')
                SdMembre = spalloc(max(Dof.NumDlA),1,0);   % Vecteur second membre
            else
                SdMembre = spalloc(max(Dof.NumDlN),1,0);   % Vecteur second membre
            end
        end
        Vect = M2*(U\(L\(M1*SdMembre(1:n0))));
        MatG.source = [Eh, zeros(m,nb); zeros(m,nh), Eb]*(MatM*...
                        (SdMembre(n0+1:end)-Mat(n0+1:end,1:n0)*Vect)) ;...
        %
        S1 = (Eh*MatM(1:nh,nh+1:nh+nb))*Hb;
        InvS2 = inv((Eb*MatM(nh+1:nh+nb,nh+1:nh+nb))*Hb);
        S = [eye(m), -S1*InvS2 ; Zero, InvS2];
        
        MatS{1} = S*[Zero,    (Eh*MatM(1:nh,1:nh))*Hh; ...
                    -eye(m),  (Eb*MatM(nh+1:nh+nb,1:nh))*Hh] ;
                
        MatS{2} = S*MatG.source;
        
%         MatS{1} = [eye(m),   +Eh*MatM(1:nh,nh+1:nh+nb)*Hb;...
%                     Zero,    +Eb*MatM(nh+1:nh+nb,nh+1:nh+nb)*Hb]\...
%                    [Zero,    +Eh*MatM(1:nh,1:nh)*Hh; ...
%                    -eye(m),  +Eb*MatM(nh+1:nh+nb,1:nh)*Hh] ;
        %
%         MatS{2} = [eye(m),   +Eh*MatM(1:nh,nh+1:nh+nb)*Hb;...
%                    Zero,    +Eb*MatM(nh+1:nh+nb,nh+1:nh+nb)*Hb]\MatG.source;
        %       
        Sdm = struct('MatM',MatM,'Eh',Eh,'Eb',Eb,'Hh',Hh,'Hb',Hb);
       
    end

end

end
%%

function [P,Q,V,Vp,InvEpz,InvMuz] = CalculVVp(Data,Phys,BetaX,BetaY,Dx,Dy)

%
% Calcul de valeurs et vecteurs propres de la matrice de transmission
%
Idx = speye(length(BetaX)); %Idx = eye(length(BetaX));
Idy = speye(length(BetaY)); %Idy = eye(length(BetaY));
%
MatP = deal(Data.nc.^2); 
MatQ = deal(ones(size(MatP)));

%
if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY)
    [MatPx,MatPy,MatPz] = deal(Data.MatP{1},Data.MatP{2},Data.MatP{3}); 
    [MatQx,MatQy,MatQz] = deal(Data.MatQ{1},Data.MatQ{2},Data.MatQ{3});
else
    [MatPx,MatPy,MatPz] = deal(Data.nc.^2); % A refaire en cas anisotropie
    [MatQx,MatQy,MatQz] = deal(ones(size(MatPx)));
end



if numel(Data.nc) == 1 && ~isfield(Data,'lx')% milieu homogène isotrope
    %
    InvEpz = 1/(Phys.K0*MatPz);
    InvMuz = 1/(Phys.K0*MatQz);
    %
    Dim = length(BetaX);
    Kxx = spdiags(BetaX(:).^2,0,Dim,Dim); %Kxx = diag(BetaX.^2);
    Kyy = spdiags(BetaY(:).^2,0,Dim,Dim); %Kyy = diag(BetaY.^2);
    Kxy = spdiags(BetaX(:).*BetaY(:),0,Dim,Dim); %Kxy = diag(BetaX.*BetaY);
    %
%     A = [1/(Phys.K0*MatP)*Kxy , Phys.K0*MatQ*Idx-1/(Phys.K0*MatP)*Kxx;...
%         -Phys.K0*MatQ*Idy+1/(Phys.K0*MatP)*Kyy ,  -1/(Phys.K0*MatP)*Kxy];

    B = [1/(Phys.K0*MatQ)*Kxy , Phys.K0*MatP*Idx-1/(Phys.K0*MatQ)*Kxx;...
        -Phys.K0*MatP*Idy+1/(Phys.K0*MatQ)*Kyy , -1/(Phys.K0*MatQ)*Kxy];

    %
    VectD = -Phys.K0^2*MatP*MatQ+(BetaX.^2)+(BetaY.^2);
    VectD = [VectD(:);VectD(:)];
    Vp = sqrt(VectD);         

    
    f = ((imag(Vp)-real(Vp))<0);
    Vp(f) = -Vp(f);
    
    %
    Ate = sqrt(-2*1i*Phys.K0*MatQx./(Dx*Dy*Vp(1:length(Vp)/2).^3));
    Atm = sqrt(2*1i./(Phys.K0*MatPx*Dx*Dy*Vp(1:length(Vp)/2))); 
    % Normalisation
    ModBeta = sqrt(BetaX.^2+BetaY.^2);
    [nx,ny] = deal(BetaX./ModBeta,BetaY./ModBeta);
    nx(isnan(nx)) = cos(Data.Phi0);
    ny(isnan(ny)) = sin(Data.Phi0);
% attention modif    
    nx(isnan(nx)) = cos(Data.Phi0);
    ny(isnan(ny)) = sin(Data.Phi0);
    %
    [Vxte,Vyte] = deal(-Ate.*ny(:),Ate.*nx(:));
    [Vxtm,Vytm] = deal(Atm.*nx(:),Atm.*ny(:));

    %V = [diag(Vxte) diag(Vxtm); diag(Vyte) diag(Vytm)];
    Dim = length(Vp);
    V = spdiags([Vxte(:); Vytm(:)],0,Dim,Dim); 
    V = spdiags([Vxtm(:);Vxtm(:)],Dim/2,V); 
    V = spdiags([Vyte(:);Vyte(:)],-Dim/2,V); 
    %
    P = B*V; %(sparse(B)*V);
else
    InvEpz = (1/Phys.K0*TfFct3D(Data,1./MatPz,0,0));
    Muy = full(Phys.K0*TfFct3D(Data,MatQy,1,0));
    Mux = full(Phys.K0*TfFct3D(Data,MatQx,0,1));
    InvMuz = full(1/Phys.K0*TfFct3D(Data,1./MatQz,0,0));
    Epy = (Phys.K0*TfFct3D(Data,MatPy,1,0));
    Epx = (Phys.K0*TfFct3D(Data,MatPx,0,1));
   
%     Kx = diag(BetaX);
%     Ky = diag(BetaY);
%     A = [Kx*InvEpz*Ky , Muy-Kx*InvEpz*Kx ; -Mux+Ky*InvEpz*Ky , -Ky*InvEpz*Kx];
%     B = [Kx*InvMuz*Ky , Epy-Kx*InvMuz*Kx ; -Epx+Ky*InvMuz*Ky , -Ky*InvMuz*Kx];
    %
    Kx1 = repmat(BetaX(:),1,length(BetaX(:)));
    Ky1 = repmat(BetaY(:),1,length(BetaY(:)));
    Kx2 = repmat(reshape(BetaX(:),1,length(BetaX(:))),length(BetaX(:)),1);
    Ky2 = repmat(reshape(BetaY(:),1,length(BetaY(:))),length(BetaY(:)),1);
    Dim = length(Kx1);

    InvEpzKx2 = InvEpz.*Kx2; InvEpzKy2 = InvEpz.*Ky2;
    %A = [Kx1.*(InvEpzKy2) , Muy-Kx1.*(InvEpzKx2) ; -Mux+Ky1.*(InvEpzKy2) , -Ky1.*(InvEpzKx2)];
    A = blkdiag(Kx1.*(InvEpzKy2) ,  -Ky1.*(InvEpzKx2)); 
    A(1:Dim,Dim+1:end) = Muy-Kx1.*(InvEpzKx2);
    A(Dim+1:end,1:Dim) = -Mux+Ky1.*(InvEpzKy2);
    %
    InvMuzKx2 = InvMuz.*Kx2; InvMuzKy2 = InvMuz.*Ky2;
    %B = [Kx1.*(InvMuzKy2) , Epy-Kx1.*(InvMuzKx2) ; -Epx+Ky1.*(InvMuzKy2) , -Ky1.*(InvMuzKx2)];
    B = blkdiag(Kx1.*(InvMuzKy2) , -Ky1.*(InvMuzKx2));
    B(1:Dim,Dim+1:end) = Epy-Kx1.*(InvMuzKx2);
    B(Dim+1:end,1:Dim) = -Epx+Ky1.*(InvMuzKy2);
    %M = [zeros(size(A,1),size(B,2)) , A ; B , zeros(size(B,1),size(A,2))];
    %clear Epx Epy Mux Muy Kx1 Kx2 Ky1 Ky2 InvEpzKx2 InvEpzKy2 InvMuzKx2 InvMuzKy2
    %
    if sum(Data.Sym) ~= 4
        [R_E,R_H,InvR_E,InvR_H] = MatSym(Data);
        A = InvR_E*A*R_H;
        B = InvR_H*B*R_E;
    end
    if isfield(Data,'Nsub') && ~isempty(Data.Nsub) && Data.Nsub ~= 0
        [P,Q,V,Vp] = deal(A,B,Data,Phys);
        return
    end

    M = A * B;
    [V,D] = eig(M,'nobalance');
    Vp = sqrt(diag(D));
    f = ((imag(Vp)-real(Vp))<0);
    Vp(f) = -Vp(f);
    P = B*V;
end    
%
%Q = V*diag(Vp);
MatVp = repmat(reshape(Vp,1,length(Vp)),length(Vp),1);
Q = V.*MatVp;

end

%%
function [R_E,R_H,InvR_E,InvR_H] = MatSym(Data)

% Consideration of symmetries in y

% Date of the latest version : 05 May 2023 (Optimized mesh)
% Author : Mondher Besbes with the help of JP Hugonin (LCF / CNRS / IOGS)

nx = 2*Data.mx+1;
ny = 2*Data.my+1;
n = nx*ny;
%
ny1 = -Data.my;
ny2 = Data.my;
%
if Data.Sym(2) ~= 2 && Data.Sym(1) == 2
    %symmetry in y 
    r0 = sparse(n,n);
    f0 = zeros(2,n);
    k = 0;
    for i0 = 1:nx
        for jp = ceil((ny1+ny2)/2):ny2
            %
            jm = ny1+ny2-jp;
            kp = i0+nx*(jp-ny1);    
            km = i0+nx*(jm-ny1);
            %
            k=k+1;
            r0(kp,k) = 1; r0(km,k) = 1; f0(k) = 1;  
            %
            if jp ~= jm
                k = k+1;
                r0(kp,k) = 1; r0(km,k) = -1; f0(k) = -1;
            end   
        end
    end
    %
    fs = find(f0==1);
    fa = find(f0==-1);
    ns = length(fs);
    na = length(fa);
    n1 = ns+na;
    %
    r=[[r0(:,fs) ; sparse(n,ns)] , [sparse(n,na) ; r0(:,fa)],...
       [r0(:,fa) ; sparse(n,na)] , [sparse(n,ns) ; r0(:,fs)]];
    
    r1=[[r0(:,fa) ; sparse(n,na)] , [sparse(n,ns) ; r0(:,fs)],...
        [r0(:,fs) ; sparse(n,ns)] , [sparse(n,na) ; r0(:,fa)]];
    %
    rr=r.';
    NormR = 1./sum(abs(rr).^2,1); 
    L = length(NormR);
    rr = rr*spdiags(NormR(:),0,L,L);
    %
    rr1=r1.';
    NormR = 1./sum(abs(rr1),1); 
    L = length(NormR);
    rr1 = rr1*spdiags(NormR(:),0,L,L);
    %
    switch Data.ChampInc
        case 1
            if Data.Sym(2)==0, f = 1:n1; end      % TM
            if Data.Sym(2)==1, f = n1+1:2*n; end  % TE
            R_E = r(:,f);
            InvR_E = rr(f,:);
            R_H = r1(:,f);
            InvR_H = rr1(f,:);
        case -1
            if Data.Sym(2)==1, f = 1:n1; end      % TE
            if Data.Sym(2)==0, f = n1+1:2*n; end  % TM
            R_H = r(:,f);
            InvR_H = rr(f,:);
            R_E = r1(:,f);
            InvR_E = rr1(f,:);
    end
else
    error('Symmetry in y is only considered, SymY=0:TM or =1:TE')
end
%
end


%%
function Tf = TfFct3D(Data,u,cx,cy)

% Calculation of matrix Tf u (permittivity / permeability)
%       tf(g)=Tf*tf(f)  ....    g=u*f
% cx: f continuous in x ....  cy: f continuous in y

% Date of the latest version : 05 May 2023 (Optimized mesh)
% Author : Mondher Besbes with the help of JP Hugonin (LCF / CNRS / IOGS)

nx = 2*Data.mx+1;
ny = 2*Data.my+1;
n = nx*ny;
%
% case of homogeneous medium
if length(unique(u(:))) == 1; Tf = u(1)*speye(n); return; end 
% 
% case of non-homogeneous medium
%
if Data.mx == 0, sens = 1; else, sens = 0; end
%
x = Data.lx(:)'; 
y = Data.ly(:)'; 
%
mx = size(x,2);
my = size(y,2);	
%
alx = (-nx+1:nx-1)*(2*pi);
aly = (-ny+1:ny-1)*(2*pi);
%
if sens==0, sens = cx|(~cy); end  % permet d'eviter l' inversion finale (sauf si cx=cy=0)
%
topx = toeplitz(nx:2*nx-1,nx:-1:1);
topy = toeplitz(ny:2*ny-1,ny:-1:1);
%

if sens == 1 % de X en Y
    if cy
        f1 = Tfu(y,u,aly);
    else
        f1 = Tfu(y,1./u,aly);
    end

    f2 = zeros(ny,ny,mx);
    if (cx&&cy) || ((~cx)&&(~cy)) 
        for ii = 1:mx, f0 = f1(ii,:); f2(:,:,ii) = f0(topy); end 
    else
        for ii = 1:mx, f0 = f1(ii,:); f2(:,:,ii) = inv(f0(topy)); end
    end

    f3 = reshape(Tfu(x,reshape(f2,ny^2,mx),alx),ny,ny,2*nx-1);
    f3 = reshape(f3(:,:,topx),[size(f3,1),size(f3,2),size(topx)]);
    f4 = reshape(permute(f3,[3, 1, 4, 2]),n,n);

    if ~cx, f4 = inv(f4); end
    
else  % de Y en X
    if cx
        f1 = Tfu(x,u.',alx);
    else
        f1 = Tfu(x,1./u.',alx);
    end
    %
    f2 = zeros(nx,nx,my);
    if (cx&&cy) || ((~cx)&&(~cy))
        for ii = 1:my, f0 = f1(ii,:); f2(:,:,ii) = f0(topx); end
    else
        for ii = 1:my, f0 = f1(ii,:); f2(:,:,ii) = inv(f0(topx)); end
    end
    f3 = reshape(Tfu(y,reshape(f2,nx^2,my),aly),nx,nx,2*ny-1);
    f3 = reshape(f3(:,:,topy),[size(f3,1),size(f3,2),size(topy)]);
    f4 = reshape(permute(f3,[1, 3, 2, 4]),n,n);
    %
    if ~cy, f4 = inv(f4); end 
      
end          

Tf = f4; %sparse(f4);

end


%-------------------------------------------------------------------------%

function tf = Tfu(x,f,alpha)

%  Fourier transform of f
%  tf = int(f*exp(-1i*alpha*x), x from 0 to 1

x = x(:)';
alpha = alpha(:).';        
tf = ([f(:,2:end),f(:,1)]-f)*exp(-1i*x(:)*alpha);
%
% % division by 1i*alpha
P0 = alpha==0;
tf(:,~P0) = -1i*tf(:,~P0)*diag(1./alpha(~P0)); 
tf(:,P0) = f*(x-[x(end)-1,x(1:end-1)]).'; 

end



