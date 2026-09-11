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

if nargin <= 3, NumLayer = []; end % FMM calculation

if isempty(NumLayer)
    if nargin == 2
        Sh = MatS;
        CoefD = CalculCoefD(Sb,[],Sh,NumLayer);
    else
        CoefD = CalculCoefD(Sb,MatS,Sh,NumLayer);
    end
else
    Sb1 = Sb; for k = 1:NumLayer(1)-1, Sb1 = ProdMatS(Sb1,MatS(k,:)); end
    if NumLayer(end) < size(MatS,1)
        Sh1 = MatS(NumLayer(end)+1,:); 
        for k = NumLayer(end)+2:size(MatS,1), Sh1 = ProdMatS(Sh1,MatS(k,:)); end
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
    if isempty(MatS) 
        MatD = ProdMatS(Sb,Sh);
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
    %
    Data = MatS{1,3};
    Phys = MatS{1,4};
    %
    %[Pdi,Pd,Pui,Pu,InvSb12] = deal(Sb{7},Sb{8},Sh{7},Sh{8},Sb{9});
    [Pdi,Pd,Pui,Pu] = deal(Sb{7},Sb{8},Sh{7},Sh{8});
    
    m = size(Sb{6},1); %size(Sb{3},1);
    
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
    if sum(Data.Sym) ~= 4, [Ib,Ih] = deal(zeros(2*m,1)); end

    mx = m/2;

    % Incidence
    if Data(1).ChampInc == +1
        if ~isempty(Phys(1).PmlX) || ~isempty(Phys(1).PmlY) 
            Ih(Sh{7}(1),1) = 1;
            Ih(Sh{7}(2),2) = 1;
        else
            if sum(Data.Sym) ~= 4
                Ih(Pui) = 1;
            else
            Ih((mx-1)/2+1+3*mx,1) = 1; %Pui = [(length(mx)-1)/2+1+length(mx)+n (length(mx)-1)/2+1+n];
            Ih((mx-1)/2+1+2*mx,2) = 1;
            end
        end
    else
        if ~isempty(Phys(1).PmlX) || ~isempty(Phys(1).PmlY) 
            Ib(Sb{7}(1),1) = 1; 
            Ib(Sb{7}(2),2) = 1;
        else
            if sum(Data.Sym) ~= 4
                Ib(Pdi) = 1;
            else
            Ib((mx-1)/2+1+mx,1) = 1; %Pdi = [(length(mx)-1)/2+1+length(mx) (length(mx)-1)/2+1];
            Ib((mx-1)/2+1,2) = 1;
            end
        end
        %
    end
    %
    if iscell(MatS{1,1}), n = 2*length(MatS{1,1}{1}); else, n = length(MatS{1,1});end
    if isempty(Ib(Pdi,:)), Fb = zeros(n,size(Ib,2)); else, Fb = -(Sb12\(Sb11*Ib(Pdi,:))); end
    if isempty(Ih(Pui,:)), Fh = zeros(n,size(Ih,2)); else, Fh = -(Sh22*Ih(Pui,:)); end
    
% Calcul de matrices et second membre
    if isfield(Data,'POD')
        if size(MatS,1) == 1
            %
            hc = Data.hc;
            Np = Data.Nsub;
            %hc = max(Mesh.CoorN(:,3))-min(Mesh.CoorN(:,3));
            dz = hc/Np;
            %
            if iscell(MatS{1,1})
                M = ProdMatCell(MatS{1,1},MatS{1,2});
                M = [[M{1,1}],[M{1,2}];[M{2,1}],[M{2,2}]];
                MatS{1,1} = cell2mat(MatS{1,1}); 
            else
                M = MatS{1,1}*MatS{1,2};
            end
            %
            n = length(MatS{1,1}); %length(A); 
            Id = speye(n);
            %
            Fb0 = dz*(MatS{1,1}*Fb); %dz*(A*Fb);
            Fh0 = dz*(MatS{1,1}*Fh); %dz*(A*Fh);
            %
            %
            if iscell(MatS{1,1}), A = cell2mat(MatS{1,1}); else, A = MatS{1,1}; end
            Ptem = Data.POD; % V*S;
            %
            for kp = 1:2
                %
            Pb = Ptem{kp}(1:n,:);
            Ph = Ptem{kp}((1:n)+Np*n,:);
            P = cell(Np,1); for k = 1:Np-1, P{k} = Ptem{kp}((1:n)+k*n,:); end
            P{Np} = Ph;
            %
            if Np ==1
                %Mtem = Pb'*(-Cb*Pb+Ph) + Ph'*(Pb-Ch*Ph);
                Mtem = Pb'*(-(Pb+dz^2/2*(M*Pb)+dz*A*(Sb12\Pb))+Ph)+...
                       Ph'*((Pb-Ph)-(dz^2/2*(M*Ph)-dz*A*(Sh21*Ph)));
            else
                %Mtem = Pb'*(-Cb*Pb+P{1}) + P{1}'*(Pb-C*P{1}+P{2});
                Mtem = Pb'*(-(Pb+dz^2/2*(M*Pb)+dz*A*(Sb12\Pb))+P{1})+...
                       P{1}'*((Pb+P{2}-2*P{1})-dz^2*(M*P{1}));
                
                for k = 2:Np-1
                    %Mtem = Mtem + P{k}'*(P{k-1}-C*P{k}+P{k+1});
                    Mtem = Mtem + P{k}'*((P{k-1}+P{k+1}-2*P{k})-dz^2*(M*P{k})); 
                end
                %Mtem = Mtem + Ph'*(P{Np-1}-Ch*Ph);
                Mtem = Mtem + Ph'*((P{Np-1}-Ph)-(dz^2/2*(M*Ph)-dz*A*(Sh21*Ph)));
    
            end
            if kp == 1
                Mtm = Mtem; 
                Ftm = Pb'*Fb0(:,1) + Ph'*Fh0(:,1);
                Etm = pinv(Mtm)*Ftm; %Mtm\Ftm;
    
            else 
                Mte = Mtem; 
                Fte = Pb'*Fb0(:,2) + Ph'*Fh0(:,2);
                Ete = pinv(Mte)*Fte; %Mte\Fte;
            end
            end
            %
            E = [Ptem{1}*Etm Ptem{2}*Ete];
            %
            % Calcul CoefD ans H
            H = zeros(size(E));
            if isempty(Ib(Pdi)) 
                Db = Sb22*(Sb12\E(1:n,:)); 
                H(1:n,:) = Sb12\E(1:n,:);
            else
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
            InvA = inv(MatS{1,1});
            for k = 1:Np-1
                H((1:n)+n*k,:) = InvA*(E((1:n)+n*(k+1),:)-E((1:n)+n*(k-1),:))/(2*dz);
            end

        else
            error('POD of multilayers: under construction! ')
        end
    else
        for k = 1:size(MatS,1)
            if ~isfield(MatS{k,3},'lx') 
                MatS{k,3}.lx = []; MatS{k,3}.ly = []; 
            end
        end
        %
        [E,H,Db,Dh] = SolveTriDiag(Sb,MatS,Sh,Fb,Fh,Ib,Ih);
    end
    CoefD = full([Dh;Db]);
end
%
end



