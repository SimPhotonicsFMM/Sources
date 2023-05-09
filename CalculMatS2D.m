function [MatS,Sdm] = CalculMatS2D(Data,Mesh,Phys,Dof,Mat,SdMembre,Sdm0)

% CalculMatS2D 
%   Calculation of S-Matrix MMF or HYB model for 2D Geometry (1D Grating)
%
% Syntax
%   MatS = CalculMatS2D(Data,Mesh,Phys);
%   MatS = CalculMatS2D(Data,Mesh,Phys,Coef);
%   MatS = CalculMatS2D(Data,Mesh,Phys,Dof,Mat);
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
if isfield(Data,'mx') && isfinite(Data.mx), mx = Data.mx; else mx = 0; end

pol = Phys.TypePol;
m = 2*mx+1;
Dx = max(Mesh.CoorN(:,1))-min(Mesh.CoorN(:,1));

Beta0 = Phys.Kx;
BetaX = Beta0 + (-mx:1:mx)*2*pi/Dx;

%% Calcul des matrices S

if Coef == -1
% --------------- côté bas -----------------------------------------------
    if ~isempty(Phys.PmlX)
        MatS = CalculMatS2D(Data,Mesh,Phys);
        P = MatS{3}; InvP = inv(P); 
        Q = MatS{4}; %InvQ = inv(Q);
        MatS{1} = [2*Q -Q*InvP; eye(m) -InvP];
        MatS{2} = zeros(length(P)*2,1);
        Vp = MatS{6};
        Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
        Pd = Pvp.';
%        Pd = 1:m;
        Pdi = [];
        MatS{1} = MatS{1}([1:m,Pd+m],[Pdi,m+1:end]);
        MatS{2} = MatS{2}([1:m,Pd+m]);
        [MatS{7},MatS{8}] = deal(Pdi,Pd);
        return
    end
    %
    Data.nc = Data.nb;
    Data.Alfa = [];
    %
    [P,Q,V,Vp] = CalculVVp(Data,Phys,BetaX,Dx,mx);
    InvP = inv(P);
    InvQ = inv(Q);
    %
    if Data.Sym(1) ~= 2
        [R,InvR] = MatSym(Data,BetaX);
        norm = sqrt(sum(R.^2,1));
        P = diag(1./norm)*InvR*P*R;
        InvP = diag(norm)*InvR*InvP*R;
        Q = diag(1./norm)*InvR*Q*R;
        InvQ = diag(norm)*InvR*InvQ*R;
        Vp = diag(InvR*diag(Vp)*R);
    end

    %
    MatSb = {[2*Q , -Q*InvP ; eye(size(Q)) , -InvP], zeros(length(P)*2,1)};
    InvSb12 = -P*InvQ;
    %
    if isfield(Data,'Eps')
        Pvp = find(abs(real(Vp))<Data.Eps); % troncature aux ordres propagatifs
    else
        Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
    end
    Pd = Pvp.';
    %
%    if Data.ChampInc == -1, Pdi = mx+1; elseif Data.ChampInc == +1, Pdi = []; end
    if Data.ChampInc == -1, Pdi = mx+1; else Pdi = []; end
    %
    n = size(MatSb{1},1)/2;
    Sb{1} = MatSb{1}([1:n,Pd+n],[Pdi,n+1:end]);
    Sb{2} = MatSb{2}([1:n,Pd+n]);
    [Sb{3},Sb{4},Sb{5},Sb{6},Sb{7},Sb{8},Sb{9}] = deal(P,Q,V,Vp,Pdi,Pd,InvSb12);
    MatS = Sb;

elseif Coef == +1
% --------------- côté haut ------------------------------------------------
    if ~isempty(Phys.PmlX)
        MatS = CalculMatS2D(Data,Mesh,Phys);
        P = MatS{3}; %InvP = inv(P);
        Q = MatS{4}; InvQ = inv(Q);
        MatS{1} = [InvQ -eye(m); P*InvQ -2*P];
        MatS{2} = zeros(length(P)*2,1);
        Vp = MatS{6};
        Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
        Pu = Pvp.';
%        Pu = 1:m;
        Pui = [];
        MatS{1} = MatS{1}([Pu,m+1:2*m],[1:m,Pui]);
        MatS{2} = MatS{2}([Pu,m+1:2*m]);
        [MatS{7},MatS{8}] = deal(Pui,Pu);
        return
    end
    %
    Data.nc = Data.nh;
    Data.Alfa = [];

    [P,Q,V,Vp] = CalculVVp(Data,Phys,BetaX,Dx,mx);

    %InvP = inv(P);
    InvQ = inv(Q);
    %
    if Data.Sym(1) ~= 2
        [R,InvR] = MatSym(Data,BetaX);
        norm = sqrt(sum(R.^2,1));
        P = diag(1./norm)*InvR*P*R;
        %InvP = diag(norm)*InvR*InvP*R;
        Q = diag(1./norm)*InvR*Q*R;
        InvQ = diag(norm)*InvR*InvQ*R;
        Vp = diag(InvR*diag(Vp)*R);
    end

    MatSh = {[InvQ , -eye(size(Q)); P*InvQ , -2*P], zeros(length(P)*2,1)};

    if isfield(Data,'Eps')
        Pvp = find(abs(real(Vp))<Data.Eps); % troncature aux ordres propagatifs
    else
        Pvp = find(abs(real(Vp))<1e-6); % troncature aux ordres propagatifs
    end
    Pu = Pvp.';

    n = size(MatSh{1},1)/2;
%    if Data.ChampInc == -1, Pui = []; elseif Data.ChampInc == +1, Pui = mx+1+n; end
    if Data.ChampInc == +1, Pui = mx+1+n; else Pui = []; end

    Sh{1} = MatSh{1}([Pu,n+1:2*n],[1:n,Pui]);
    Sh{2} = MatSh{2}([Pu,n+1:2*n]);
    [Sh{3},Sh{4},Sh{5},Sh{6},Sh{7},Sh{8}] = deal(P,Q,V,Vp,Pui,Pu);
    MatS = Sh;

else
% --------------- couche homogène ou réseau ------------------------------------------------


    if nargin == 3,     % MMF
                
        
        Pn = find(abs(Mesh.CoorN(:,2)-min(Mesh.CoorN(:,2)))<Dx/1e6);
        x = sort(Mesh.CoorN(Pn,1)-min(Mesh.CoorN(:,1)));
        dx = cumsum(diff(x))';
        Data.Alfa = dx(1:end-1)/Dx; % à modifier pour le calcul du champ

        Pe = unique(find(sum(ismember(Mesh.Cn,Pn),2)));
        [~,P] = sort(Mesh.CoorF(Pe,1)); 
        Data.nc = Data.Indice(Mesh.Nsd(Pe(P))); 
        li = diff(x);
        if ~isempty(Data.Alfa)
            NumSD = Mesh.Nsd(Pe(P))';
            if ~any(abs(x-Dx/2)<Dx/1e6)
                x = sort([x; Dx/2]);
                P = find(abs(x-Dx/2)<Dx/1e6);
                NumSD = [NumSD(1:P-1) NumSD(P-1) NumSD(P:end)];
                li = diff(x);
            end
            xa = x(1:end-1)+diff(x)/2;
            P=xa>Dx/2;
            Num = [NumSD(P(1:end)) NumSD(~P)];
            lx = [li(P(1:end)); li(~P)];
            ff = cumsum(lx)/Dx;
            ff = ff(1:end-1)';
            %
            Data.nc = Data.Indice(Num);
            Data.Alfa = ff;
            %disp(Data.nc)
            %disp(Data.Alfa)
        end

        [P,Q,V,Vp,TpzMatQy] = CalculVVp(Data,Phys,BetaX,Dx,mx);
        InvP = inv(P);
        InvQ = inv(Q);
        %InvP = InvMat(P);
        %InvQ = InvMat(Q);
        %
        hc = max(Mesh.CoorN(:,2))-min(Mesh.CoorN(:,2)); %Data.hc;
        %
        MatS{1} = [Q*diag(1./cosh(Vp*hc))*InvQ , Q*diag(tanh(Vp*hc))*InvP;
                  -P*diag(tanh(Vp*hc))*InvQ    , P*diag(1./cosh(Vp*hc))*InvP];
        MatS{2} = zeros(length(P)*2,1);
        [MatS{3},MatS{4},MatS{5},MatS{6},MatS{7}] = deal(P,Q,V,Vp,TpzMatQy);
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
function [P,Q,V,Vp,TpzMatQy] = CalculVVp(Data,Phys,BetaX,Dx,mx)

pol = Phys.TypePol;

if ~isfield(Data,'nc'), Data.nc = Data.Indice; end
if ~isfield(Data,'Alfa'), Data.Alfa = []; end

[MatPx,MatPy,MatPz] = deal(Data.nc.^(2-pol)); % Formulation en E (nodes)
[MatQx,MatQy,MatQz] = deal(Data.nc.^(pol));

if ~isempty(Phys.PmlX)
    Eps = Phys.CaractEps;
    Nu = Phys.CaractNu;  % en 2D EF Nu1 
    [MatPx,MatPy,MatPz] = deal(Eps(1,1,:),Eps(2,2,:),Eps(3,3,:));
    [MatQx,MatQy,MatQz] = deal(1./Nu(1,1,:),1./Nu(2,2,:),1./Nu(3,3,:));
end

%
K = 2*pi*(-2*mx:1:2*mx);
X = [Data.Alfa 1];
%
TfMatQy = TfFct(MatQy,X,K);%            
TpzMatQy = toeplitz(TfMatQy(2*mx+1:end),TfMatQy(2*mx+1:-1:1));
%
TfMatPz = TfFct(MatPz,X,K);%            
TpzMatPz = toeplitz(TfMatPz(2*mx+1:end),TfMatPz(2*mx+1:-1:1));
%
TfInvMatQx = TfFct(1./MatQx,X,K);%            norm(TfMatQ0-TfMatQ);
TpzInvMatQx = toeplitz(TfInvMatQx(2*mx+1:end),TfInvMatQx(2*mx+1:-1:1));
%
A = Phys.K0*inv(TpzInvMatQx);
B = -Phys.K0*TpzMatPz + 1/Phys.K0*diag(BetaX)*inv(TpzMatQy)*diag(BetaX);
%
if isfield(Data,'DofType') && strcmp(Data.DofType,'Edge')
    [A,B] = deal(B,A);
end




if length(Data.nc) > 1 
    % Matrice Valeurs Propres
    if Data.Sym(1) ~= 2
        [R,InvR] = MatSym(Data,BetaX);
        A = InvR*A*R;
        B = InvR*B*R;
    end
    [V,D] = eig(A*B);
    Vp = sqrt(diag(D));
    f = find((imag(Vp)-real(Vp))<0); 
    Vp(f)=-Vp(f);
else
    % Matrice Valeurs Propres
    D = diag(-Phys.K0^2*MatQx*MatPz + MatQx/MatQy*BetaX.^2);
    Vp = sqrt(diag(D));
    f = find((imag(Vp)-real(Vp))<0); 
    Vp(f)=-Vp(f);
    % Matrice Vecteurs Propres (Flux de Poyntong = 1)
    V = diag(sqrt(-2i*Phys.K0*MatQx/Dx./Vp.^3));

    if isfield(Data,'DofType') && strcmp(Data.DofType,'Edge')
        V = diag(sqrt(-2i/Phys.K0/MatQx/Dx./Vp));
    end

end

Q = V*diag(Vp);
P = B*V; 

%[ cond(TpzMatQy) cond(P) cond(Q)],

end

%%
function [R,InvR] = MatSym(Data,BetaX)

N = length(BetaX);
Ns = ceil(N/2);
Na = N-Ns;
%
R = sparse(1:N,N:-1:1,ones(1,N),N,N)+sparse(1:N,1:N,[ones(1,Ns) -ones(1,Na)],N,N);
InvR = inv(R);
%
x0 = 0;  % Géométrie symétrie par rapport à zéro
R = diag(exp(-1i*BetaX*x0))*R;
InvR = InvR*diag(exp(+1i*BetaX*x0));
%
if Data.TypePol == 0, P = 1:Ns; else, P = Ns+1:N; end
%
R = R(:,P);
InvR = InvR(P,:);

end
%%
function Tf = TfFct(Fct,X,Beta)

% TF d'une fonction constante par intervale

X = X(:);
Fct = Fct(:).';
Beta = Beta(:).';

Tf = ((Fct-[Fct(2:end),Fct(1)])*exp(-1i*X*Beta)).*(1i./Beta);

Pm = Beta == 0;
Tf(Pm) = Fct * (X-[0;X(1:end-1)]);

%if length(Beta) == 1, return, end

%Tf(~Pm) = (Fct-[Fct(2:end),Fct(1)])*exp(-1i*X*Beta(~Pm))*diag(1i./Beta(~Pm));

end
