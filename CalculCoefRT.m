function [r,t,CoefD,R,T,E,H] = CalculCoefRT(Sb,MatS,Sh,NumLayer)

% CalculCoefRT
%   Calculation of diffraction efficiency by a product of S-matrices
%
% Syntax
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh);
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,Sh);
%   [r,t,CoefD,R,T,E,H] = CalculCoefRT(Sb,MatS,Sh,NumLayer);
%
% Description
%   Sb   : Substrat S-matrix
%   MatS : S-Matrices of different layers
%   Sh   : Upstrat S-matrix
%
%   r : Reflectivity
%   t : Transmittivity
%   CoefD : Complex diffraction coeficients
%   R,T : Fresnel coefficients
%
% Example : Gold grating+substrat
%   ld = .85;
%   dx = .2; lix = .1; dy = .2; liy = .1; h = .03;
%   Mesh = MeshLayer(dx,lix,dy,liy,h,2,2,2); 
%   Data = SetData('Lambda0',ld,'Theta0',0,'Phi0',0,'ChampInc',-1,'TypePol',2,...
%       'mx',5,'my',5,'nh',1.33,'nb',1.7,'Indice',[IndexVal('Au',ld) 1.33]);
%   Phys = CaractMat(Mesh,Data);
%   Sb = CalculMatS(Data,Mesh,Phys,-1); % milieu bas
%   Sh = CalculMatS(Data,Mesh,Phys,+1); % milieu haut
%   MatS = CalculMatS(Data,Mesh,Phys); 
%   [r,t,CoefD,R,T] = CalculCoefRT(Sb,MatS,Sh);
%

% Date of the latest version : 13 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

% m = size(Sb{3},1);
% 
% if nargin == 2 
%     S = Sb;
%     Sh = MatS;
%     MatD = ProdMatS(S,Sh);
% else    
%     %S = ProdMatS(Sb,ProdMatS(MatS));
%     S = Sb; for k=1:size(MatS,1), S = ProdMatS(S,MatS(k,:)); end
%     MatD = ProdMatS(S,Sh);
% end
% 
% %
% %MatD = ProdMatS(S,Sh);
% %
% CoefD = MatD{1};
% if isempty(CoefD), CoefD = MatD{2}; end

if nargin == 3, NumLayer = []; end % FMM calculation

if isempty(NumLayer)
    CoefD = CalculCoefD(Sb,MatS,Sh,NumLayer);
else
    Sb1 = Sb; for k = 1:NumLayer-1, Sb1 = ProdMatS(Sb1,MatS(k,:)); end
    if NumLayer < size(MatS,1)
        Sh1 = MatS(NumLayer+1,:); 
        for k = NumLayer+2:size(MatS,1), Sh1 = ProdMatS(Sh1,MatS(k,:)); end
        Sh1 = ProdMatS(Sh1,Sh);
    else
        Sh1 = Sh;
    end

    [CoefD,E,H] = CalculCoefD(Sb1,MatS(NumLayer,:),Sh1,NumLayer);
end

m = size(Sb{3},1);

[Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});
%
Dim = size(CoefD,2);
%
if Dim == 2 %|| numel(MatS{1}) == numel(MatS{2})
    [r,t] = deal(zeros(m,Dim));
    [R,T] = deal(zeros(m,Dim));
    %
    rt = abs(CoefD).^2;
    
    t(Pu,1) = rt(1:length(Pu),1)';
    r(Pd,1) = rt(length(Pu)+1:end,1)';
    %
    T(Pu,1) = CoefD(1:length(Pu),1).';
    R(Pd,1) = CoefD(length(Pu)+1:end,1).';
else
    %
    rt = abs(CoefD).^2;
    
    t = sum(rt(1:length(Pu),1));
    r = sum(rt(length(Pu)+1:end,1));
    %
    T = CoefD(1:length(Pu),1);
    R = CoefD(length(Pu)+1:end,1);

end


%
if size(rt,2) == 2 
    t(Pu,2) = rt(1:length(Pu),2)'; 
    r(Pd,2) = rt(length(Pu)+1:end,2)';
    %
    T(Pu,2) = CoefD(1:length(Pu),2).'; 
    R(Pd,2) = CoefD(length(Pu)+1:end,2).';

end

if isempty(Pdi), [r,t] = deal(t,r); [R,T] = deal(T,R); end  % cas incidence en bas 

end

% ------------------------------------------------------------------------%

function [CoefD,E,H] = CalculCoefD(Sb,MatS,Sh,NumLayer)

if isempty(NumLayer)
    if nargin == 2 
        S = Sb;
        Sh = MatS;
        MatD = ProdMatS(S,Sh);
    else    
        %S = ProdMatS(Sb,ProdMatS(MatS));
        S = Sb; for k=1:size(MatS,1), S = ProdMatS(S,MatS(k,:)); end
        MatD = ProdMatS(S,Sh);
    end
    
    %
    %MatD = ProdMatS(S,Sh);
    %
    CoefD = MatD{1};
    if isempty(CoefD), CoefD = MatD{2}; end
    

else  % Calcul numérique par schéma semi-implicite

    Data = MatS{3};
    Phys = MatS{4};
    hc = Data.hc;

    A = MatS{1};
    B = MatS{2};
    %
    %[Pdi,Pd,Pui,Pu,InvSb12] = deal(Sb{7},Sb{8},Sh{7},Sh{8},Sb{9});
    [Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});
    
    m = size(Sb{3},1);
    
    Sh1 = Sh{1};
    [n1,n2,n3,n4] = deal(m,length(Pui),length(Pu),m);
    [Sh11,Sh12,Sh21,Sh22] = deal(Sh1(1:n3,1:n1),Sh1(1:n3,n1+1:end),...
                                 Sh1(n3+1:end,1:n1),Sh1(n3+1:end,n1+1:end));
    Sb1 = Sb{1};
    [n1,n2,n3,n4] = deal(length(Pdi),m,m,length(Pd));
    [Sb11,Sb12,Sb21,Sb22] = deal(Sb1(1:n3,1:n1),Sb1(1:n3,n1+1:end),...
                                 Sb1(n3+1:end,1:n1),Sb1(n3+1:end,n1+1:end));
    %
    %
    [Ib,Ih] = deal(zeros(2*m,2)); % % dim 2 A vérifier
    
    mx = m/2;

    % Incidence
    if Data.ChampInc == +1
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY) 
            Ih(Sh{7}(1),1) = 1;
            Ih(Sh{7}(2),2) = 1;
        else
            Ih((mx-1)/2+1+3*mx,1) = 1; %Pui = [(length(mx)-1)/2+1+length(mx)+n (length(mx)-1)/2+1+n];
            Ih((mx-1)/2+1+2*mx,2) = 1;
        end
    else
        if ~isempty(Phys.PmlX) || ~isempty(Phys.PmlY) 
            Ib(Sb{7}(1),1) = 1; 
            Ib(Sb{7}(2),2) = 1;
        else
            Ib((mx-1)/2+1+mx,1) = 1; %Pdi = [(length(mx)-1)/2+1+length(mx) (length(mx)-1)/2+1];
            Ib((mx-1)/2+1,2) = 1;
        end
        %
    end
% Calcul de matrices et second membre
    Np = Data.Nsub;
    %hc = max(Mesh.CoorN(:,3))-min(Mesh.CoorN(:,3));
    dz = hc/Np;
    n = length(A); %n=200;
    Id = speye(n);
    %
    %if isempty(Ib(Pdi,:)), Fb = zeros(n,2); else, Fb = -dz*ASb12*(Sb11*Ib(Pdi,:)); end
    if isempty(Ib(Pdi,:)), Fb = zeros(n,2); else, Fb = -dz*A*(Sb12\(Sb11*Ib(Pdi,:))); end
    if isempty(Ih(Pui,:)), Fh = zeros(n,2); else, Fh = -dz*A*(Sh22*Ih(Pui,:)); end
    %

    if isfield(Data,'POD')
        %
        Ptem = Data.POD; %Ptm = Data.POD{1}; % V*S;
        %Pte = Data.POD{2}; % V*S
        %
        for kp = 1:2
            %
        Pb = Ptem{kp}(1:n,:);
        Ph = Ptem{kp}((1:n)+Np*n,:);
        P = cell(Np,1); for k = 1:Np-1, P{k} = Ptem{kp}((1:n)+k*n,:); end
        P{Np} = Ph;
        %
        %Mtm = (-Pb'*Cb+P1')*Pb + (Pb'*C1-P1'*C+Ph'*C1)*P1 + (P1'-Ph'*Ch)*Ph;
        PbA = Pb'*A;
        PhA = Ph'*A;
        if Np ==1
            %Mtem = Pb'*(-Cb*Pb+Ph) + Ph'*(Pb-Ch*Ph);
            Mtem = -(Pb'*Pb+dz^2/2*(PbA)*(B*Pb)+dz*(PbA)*(Sb12\Pb))+Pb'*Ph+...
                    Ph'*(Pb-Ph)-(dz^2/2*(PhA)*(B*Ph)-dz*(PhA)*(Sh21*Ph));
        else
            %Mtem = Pb'*(-Cb*Pb+P{1}) + P{1}'*(Pb-C*P{1}+P{2});
            Mtem = -(Pb'*Pb+dz^2/2*(PbA)*(B*Pb)+dz*(PbA)*(Sb12\Pb))+Pb'*P{1}+...
                     P{1}'*(Pb+P{2}-2*P{1})-dz^2*(P{1}'*A)*(B*P{1});

            for k = 2:Np-1
                %Mtem = Mtem + P{k}'*(P{k-1}-C*P{k}+P{k+1});
                Mtem = Mtem + P{k}'*(P{k-1}+P{k+1}-2*P{k})-dz^2*(P{k}'*A)*(B*P{k}); 
            end
            %Mtem = Mtem + Ph'*(P{Np-1}-Ch*Ph);
            Mtem = Mtem + Ph'*(P{Np-1}-Ph)-(dz^2/2*(PhA)*(B*Ph)-dz*(PhA)*(Sh21*Ph));

        end
        if kp == 1
            Mtm = Mtem; 
            Ftm = Pb'*Fb(:,1) + Ph'*Fh(:,1);
            Etm = Mtm\Ftm;

        else 
            Mte = Mtem; 
            Fte = Pb'*Fb(:,2) + Ph'*Fh(:,2);
            Ete = Mte\Fte;
        end
        end
        %
        E = [Ptem{1}*Etm Ptem{2}*Ete];
        
    else
    InvSb12 = Sb12\Id; %inv(Sb12);

    ASb12 = A*InvSb12;  
    M = A*B;
    %clear B
    %
    C = 2*Id + dz^2*M;                 
    C1 = Id;% - dz^2/4*M;                
    Cb = Id + dz^2/2*M + dz*ASb12; 
    Ch = Id + dz^2/2*M - dz*A*Sh21;    

    % 
    %clear M ASb12
          
    test = 0;
    if test == 2
% Décomposition LU
    %bk = cell(Np+1,1);
    ck = cell(Np+1,1);
    %for k = 1:Np+1, ck{k} = nan(size(C)); end

    InvCb = Cb\Id;%inv(Cb);%ck{1} = -Cb;%bk{1} = -Cb;

    ck{1} = -InvCb*C1;%ck{1} = inv(ck{1})*C1;%ck{1} = inv(bk{1})*C1;%bk{1}\C1;
    %bk{1} = [];
    for k = 2:Np
        ck{k} = -C-ck{k-1};%bk{k} = -C-ck{k-1};
        ck{k} =  inv(ck{k});%ck{k} =  inv(bk{k});%bk{k}\Id;
        %bk{k} = [];
    end
    ck{Np+1} = -Ch-C1*ck{Np};%bk{Np+1} = -Ch-C1*ck{Np}; 
    %ck{Np+1} = inv(ck{Np+1});%ck{Np+1} = inv(bk{Np+1});%bk{Np+1}\Id;%
    clear Cb Ch C
% Résolution
    [xk,yk] = deal(cell(Np+1,1));
    yk{1} = -InvCb*Fb;%yk{1} = -Cb\Fb;%bk{1}\Fb;
    clear InvCb
    for k = 2:Np, yk{k} = -ck{k}*yk{k-1}; end
%    yk{Np+1} = ck{Np+1}*(Fh-C1*yk{Np});
    yk{Np+1} = ck{Np+1}\(Fh-C1*yk{Np});
    clear C1
    %
    ck{Np+1} = [];
    xk{Np+1} = yk{Np+1};
    for k = Np:-1:1, xk{k} = yk{k}-ck{k}*xk{k+1}; ck{k}=[];end
    E = nan(n*(Np+1),2);
    for k = 1:Np+1, E((1:n)+n*(k-1),:) = xk{k}; end

    else

      if Np == 1
        Mat = blkdiag(-Cb,-Ch); Mat(n+1:end,1:n) = C1; Mat(1:n,n+1:end) = C1;
        E = Mat\[Fb;Fh];
      else
          if test == 0
        InvCb = inv(Cb);%Cb\Id;
        B1 = C - InvCb*C1;
        F1 = InvCb*Fb;
        %
        for k = 2:Np-1
            InvB1 = inv(B1);%B1\Id;%
            F1 = InvB1*F1;
            B1 = C - InvB1;
        end
        %
        Mat = blkdiag(-B1,-Ch); Mat(n+1:end,1:n) = C1; Mat(1:n,n+1:end) = eye(n);
        TabE = Mat\[F1;Fh];
        %no
        xk = cell(Np+1,1);
        xk{Np} = TabE(1:n,:); % E(Np-1)
        xk{Np+1} = TabE(n+1:end,:); % Eh
        for k = Np-1:-1:1, xk{k} = C*xk{k+1} - xk{k+2}; end
            %
        E = nan(n*(Np+1),2);
        for k = 1:Np+1, E((1:n)+n*(k-1),:) = xk{k}; end
        %
          else
        TestInv;
          end
      end
    end

    clear ck
    end
% Calcul CoefD
    
    H = zeros(size(E));
    if isempty(Ib(Pdi)) 
        Db = Sb22*(Sb12\E(1:n,:)); %Sb22*InvSb12*E(1:n,:);
        H(1:n,:) = Sb12\E(1:n,:);%InvSb12*E(1:n,:);
    else
        %Db = (Sb21-Sb22*InvSb12*Sb11)*Ib(Pdi,:)+Sb22*InvSb12*E(1:n,:);
        %H(1:n,:) = InvSb12*E(1:n,:) - InvSb12*(Sb11*Ib(Pdi,:));
        Db = Sb21*Ib(Pdi,:)-Sb22*(Sb12\(Sb11*Ib(Pdi,:)))+Sb22*(Sb12\E(1:n,:));
        H(1:n,:) = Sb12\E(1:n,:) - Sb12\(Sb11*Ib(Pdi,:));
    end
    %
    
    if isempty(Ih(Pui))
        Dh = Sh11*E(n*Np+1:n*(Np+1),:);
        H(n*Np+1:n*(Np+1),:) = Sh21*E(n*Np+1:n*(Np+1),:);
    else
        Dh = Sh11*E(n*Np+1:n*(Np+1),:) + Sh12*Ih(Pui,:);
        H(n*Np+1:n*(Np+1),:) = Sh21*E(n*Np+1:n*(Np+1),:)+ Sh22*Ih(Pui,:);
    end
    %
    if ~isfield(Data,'POD')
    InvA = inv(A);
    for k = 1:Np-1
        H((1:n)+n*k,:) = (InvA*(E((1:n)+n*(k+1),:)-E((1:n)+n*(k-1),:)))/(2*dz);
    end
    end
    
    CoefD = full([Dh;Db]);
end

end

%-------------------------------------------------------------%
function [Mtm,Mte,Ftm,Fte] = CalculPOD(Ptm,Pte,C,C1,Cb,Ch,Fb,Fh,Np)

n = length(C);


Pb = Ptm(1:n,:);
P1 = Ptm((1:n)+(Np-1)*n,:);
Ph = Ptm((1:n)+Np*n,:);

Mtm = (-Pb'*Cb+P1')*Pb + (Pb'*C1-P1'*C+Ph'*C1)*P1 + (P1'-Ph'*Ch)*Ph;
Ftm = Pb'*Fb(:,1) + Ph'*Fh(:,1);

Pb = Pte(1:n,:);
P1 = Pte((1:n)+(Np-1)*n,:);
Ph = Pte((1:n)+Np*n,:);

Mte = (-Pb'*Cb+P1')*Pb + (Pb'*C1-P1'*C+Ph'*C1)*P1 + (P1'-Ph'*Ch)*Ph;
Fte = Pb'*Fb(:,2) + Ph'*Fh(:,2);

return

Mat = sparse(n*(Np+1),n*(Np+1));
[i,j] = find(Cb+eps);
Mat = Mat + sparse(i,j,-Cb,n*(Np+1),n*(Np+1));
[i,j] = find(C+eps);
for k = 2:Np 
    Mat = Mat + sparse(i+n*(k-1),j+n*(k-1),-C,n*(Np+1),n*(Np+1)); 
end
[i,j] = find(Ch+eps);
Mat = Mat + sparse(i+n*Np,j+n*Np,-Ch,n*(Np+1),n*(Np+1));
%
[i,j] = find(C1+eps);
Mat = Mat + sparse(i,j+n,C1,n*(Np+1),n*(Np+1));
Mat = Mat + sparse(i+n*Np,j+n*Np-n,C1,n*(Np+1),n*(Np+1));
Mat = Mat + sparse(n+1:n*Np,2*n+1:n*(Np+1),diag(speye(n*(Np-1))),n*(Np+1),n*(Np+1));
Mat = Mat + sparse(n+1:n*Np,1:n*(Np-1),diag(speye(n*(Np-1))),n*(Np+1),n*(Np+1));
%
F = sparse(n*(Np+1),2);
F(1:n,:)= Fb;
F(n*Np+1:n*(Np+1),:) = Fh;

end

