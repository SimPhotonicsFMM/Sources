classdef utilMesh

% utilMesh
%   class of methods to modify a mesh
%
% Properties
%   Mesh : Mesh of a structure
%
% Methods
%   AssembMesh : Add two meshes
%       Mesh = AssembMesh(Mesh1,Mesh2)
%
%   ExtraxtMesh : Extract part of a mesh specified by a list of elements Pe
%       Mesh1 = ExtractMesh(Mesh,Pe)
%
%   MoveMesh    : Move a mesh at a distance "Dep"
%       Mesh1 = MoveMesh(Mesh,Dep)
%
%   ShifMesh    : Shift a mesh for field calculation
%       Mesh1 = ShiftMesh(Mesh)
%
% Example
%   Mesh = MeshLayer(1,.5,1,.5,.2,2,2,2); % cell with a square inclusion
%   Pe = Mesh.CoorV(:,1)>0 & Mesh.CoorV(:,2)>0; % x>0 and y>0
%   Mesh1 = ExtractMesh(utilMesh(Mesh),Pe);
%   figure,VisuMesh(Mesh1)

% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

    properties
        mesh = struct();
    end
    methods
        function obj = utilMesh(val)
            if nargin == 1
                obj.mesh = val;
            end
        end
        %
        % Assemb two meshes
        function Mesh = AssembMesh(obj1,obj2)
            Mesh1 = [obj1.mesh];
            Mesh2 = [obj2.mesh];
            %
            N = size(Mesh1.CoorN,1);
            A = size(Mesh1.CoorA,1);
            F = size(Mesh1.CoorF,1);
            %
            Mesh.Cn = [Mesh1.Cn; Mesh2.Cn+N];
            Mesh.Ca = [Mesh1.Ca; abs(Mesh2.Ca)+A];
            Mesh.Cf = [Mesh1.Cf; Mesh2.Cf+F];
            Mesh.CoorN = [Mesh1.CoorN; Mesh2.CoorN];
            Mesh.CoorA = [Mesh1.CoorA; Mesh2.CoorA];
            Mesh.CoorF = [Mesh1.CoorF; Mesh2.CoorF];
            Mesh.CoorV = [Mesh1.CoorV; Mesh2.CoorV];
            Mesh.Nsd = [Mesh1.Nsd; Mesh2.Nsd];
            Mesh.ExtAr = [Mesh1.ExtAr; Mesh2.ExtAr+N];
            Mesh.ExtFa = [Mesh1.ExtFa; Mesh2.ExtFa+N];
            Mesh.ExtFaAr = [Mesh1.ExtFaAr; abs(Mesh2.ExtFaAr)+A];
            Mesh.NsdF = [Mesh1.NsdF; Mesh2.NsdF];
            Mesh.TabP = [Mesh1.TabP; Mesh2.TabP];
            Mesh.Vol = [Mesh1.Vol; Mesh2.Vol];
            Mesh.TypeElmt = Mesh1.TypeElmt;
            %
            if isfield(Mesh1,'xv') && ~isempty(Mesh1.xv)
                Mesh.xv = Mesh1.xv;
                Mesh.yv = Mesh1.yv;
            end
        end
        %
        % Extract a part of mesh
        function Mesh1 = ExtractMesh(obj,Pe)
            %
            Mesh = [obj.mesh];
            Mesh1 = Mesh;
            %
            Mesh1.Nsd = Mesh.Nsd(Pe);
            NumSD = unique(Mesh1.Nsd);
%             for k = 1:length(NumSD)
%                 Pe1 = Mesh1.Nsd == NumSD(k);
%                 Mesh1.Nsd(Pe1) = k;
%             end
            %Mesh1.Nsd = Mesh1.Nsd-min(Mesh1.Nsd)+1;
%            NumSD = unique(Mesh1.Nsd);
            %if size(Mesh.CoorN,2) == 2, Mesh1.Surf = Mesh.Surf(Pe); end
            %
            Pn = unique(Mesh.Cn(Pe,:));
            Pa = unique(abs(Mesh.Ca(Pe,:)));
            %
            Cn = Mesh.Cn(Pe,:);
            Ca = Mesh.Ca(Pe,:);
            %
            Mesh1.CoorN =  Mesh.CoorN(Pn,:);
            P1 = zeros(1,size(Mesh.CoorN,1));
            for k = 1:length(Pn), P1(Pn(k)) = k; end  % compresser la numérotation des noeuds
            Mesh1.Cn = P1(Cn);
            ExtAr = Mesh.ExtAr(Pa,:);
            Mesh1.ExtAr = P1(ExtAr);
            
            %
            if size(Mesh.CoorN,2) == 3 
                Pf = unique(Mesh.Cf(Pe,:));
                Cf = Mesh.Cf(Pe,:);
                ExtFa = Mesh.ExtFa(Pf,:);
                Mesh1.ExtFa = P1(ExtFa);
            end
            %
            Mesh1.CoorA =  Mesh.CoorA(Pa,:);
            P1 = int32(zeros(1,size(Mesh.CoorA,1)));
            for k = 1:length(Pa), P1(Pa(k)) = k; end  % compresser la numérotation des arêtes
            Mesh1.Ca = P1(abs(Ca)).*sign(Ca);
            
            if size(Mesh.CoorN,2) == 3
                ExtFaAr = Mesh.ExtFaAr(Pf,:);
                Mesh1.ExtFaAr = P1(abs(ExtFaAr)).*sign(ExtFaAr);
                %
                Mesh1.CoorF =  Mesh.CoorF(Pf,:);
                P1 = zeros(1,size(Mesh.CoorF,1));
                for k = 1:length(Pf), P1(Pf(k)) = k; end  % compresser la numérotation des arêtes
                Mesh1.Cf = P1(Cf);
                %
                Mesh1.CoorV=  Mesh.CoorV(Pe,:);
                Mesh1.NsdF = Mesh.NsdF(Pf);
                Mesh1.Vol = Mesh.Vol(Pe);
                %Mesh1.Surf = Mesh.Surf(Pf);
                Mesh1.TabP = Mesh.TabP(Pf);
                %Mesh1.VectNorm = Mesh.VectNorm(Pf,:);
                Mesh1.TabNsdF = {};
                for k = 1:length(NumSD) 
                    %Tab = int8(zeros(size(Mesh.TabNsdF{NumSD(k)})));
                    Tab = Mesh.TabNsdF{1,NumSD(k)};
                    Mesh1.TabNsdF{1,k} = int8(zeros(size(Pf)));
                    P1 = Tab(Pf) ~= 0;
                    Mesh1.TabNsdF{1,k}(P1) = NumSD(k);%Tab(Pf)+k;
                    %Mesh1.TabNsdF{1,k} = Tab(Pf)+NumSD(k);%Tab(Pf)+k; 
                end 
            end
            %
                
        end
        %
        % Move a Mesh (translation)
        function Mesh2 = MoveMesh(obj,Dep)
    
            % MoveMesh 
            %   Déplacer un maillage 2D ou 3D de Dep 
            %
            % Syntax
            %   Mesh2 = MoveMesh(Mesh1,Dep);
            %
            % Description
            %   Mesh1 : structure données du maillage initial 
            %   Dep : vecteur de translation
            %
            %   Mesh2 : structure données du maillage final après déplacement
            
            Mesh1 = [obj.mesh];
            Mesh2 = Mesh1;
            
            for kc = 1:length(Mesh1)
                % Nodes
                TabDep = repmat(Dep,size(Mesh1(kc).CoorN,1),1);
                Mesh2(kc).CoorN = Mesh1(kc).CoorN + TabDep;
                % Edges
                TabDep = repmat(Dep,size(Mesh1(kc).CoorA,1),1);
                Mesh2(kc).CoorA = Mesh1(kc).CoorA + TabDep;
                % Facets
                if isfield(Mesh1(kc),'CoorF')
                    TabDep = repmat(Dep,size(Mesh1(kc).CoorF,1),1);
                    Mesh2(kc).CoorF = Mesh1(kc).CoorF + TabDep;
                end
                %
                % Volumes
                if isfield(Mesh1(kc),'CoorV')
                    TabDep = repmat(Dep,size(Mesh1(kc).CoorV,1),1);
                    Mesh2(kc).CoorV = Mesh1(kc).CoorV + TabDep;
                end
            end
        end
        %
        % Shift a Mesh (as fftshit)
        function Mesh = ShiftMesh(obj)
    
            Mesh0 = [obj.mesh];
            [xmin,xmax,ymin,ymax] = MinMax(util(Mesh0(1).CoorN));
            Mesh = Mesh0;
            
            for k = 1:length(Mesh0)
                if max(Mesh0(k).Nsd) > 1
                    Pe = find(Mesh0(k).CoorV(:,1)>0 & Mesh0(k).CoorV(:,2)>0);
                    Mesh1 = ExtractMesh(utilMesh(Mesh0(k)),Pe);
                    Mesh1 = MoveMesh(utilMesh(Mesh1),[-xmax -ymax 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)>0 & Mesh0(k).CoorV(:,2)<0);
                    Mesh2 = ExtractMesh(utilMesh(Mesh0(k)),Pe);
                    Mesh2 = MoveMesh(utilMesh(Mesh2),[-xmax -ymin 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)<0 & Mesh0(k).CoorV(:,2)<0);
                    Mesh3 = ExtractMesh(utilMesh(Mesh0(k)),Pe);
                    Mesh3 = MoveMesh(utilMesh(Mesh3),[-xmin -ymin 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)<0 & Mesh0(k).CoorV(:,2)>0);
                    Mesh4 = ExtractMesh(utilMesh(Mesh0(k)),Pe);
                    Mesh4 = MoveMesh(utilMesh(Mesh4),[-xmin -ymax 0]);
                    %
                    if length(Mesh0) == 1
                        Mesh = AssembMesh(utilMesh(AssembMesh(utilMesh(Mesh1),utilMesh(Mesh2))),...
                                   utilMesh(AssembMesh(utilMesh(Mesh3),utilMesh(Mesh4))));
                    else
                        Mesh(k) = AssembMesh(utilMesh(AssembMesh(utilMesh(Mesh1),utilMesh(Mesh2))),...
                                   utilMesh(AssembMesh(utilMesh(Mesh3),utilMesh(Mesh4))));
                    end
                end
            end
        
        end
        % New number of subdomain
        function Mesh = SetNumSD(obj,Num0,Num1)

            Mesh0 = [obj.mesh];
            Mesh = Mesh0;
            if nargin == 2
                Pe = Num0;
                kd = max(Mesh0.Nsd)+1;
                Mesh.Nsd(Pe) = kd;
                Pf = unique(Mesh.Cf(Pe,:));
                Mesh.NsdF(Pf) = kd;
                Mesh.TabNsdF{kd}(1:length(Mesh.CoorF),1) = 0;
                Mesh.TabNsdF{kd}(Pf) = kd;
            else
            for k = 1:length(Mesh0)
                for kd = 1:length(Num0)
                    Pe = Mesh0(k).Nsd == Num0(kd);
                    Mesh(k).Nsd(Pe) = Num1(kd);
                end
                %
                %Mesh.TabNsdF = cell(max(Mesh.Nsd),1);
                for kd = Num1
                    Pe = Mesh(k).Nsd == kd;
                    Pf = unique(Mesh(k).Cf(Pe,:));
                    Mesh(k).NsdF(Pf) = kd;
                    Mesh(k).TabNsdF{kd}(1:length(Mesh(k).CoorF),1) = 0;
                    Mesh(k).TabNsdF{kd}(Pf) = kd;
                end
            end
            end
        end

    end
end

