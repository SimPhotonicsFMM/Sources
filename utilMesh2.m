classdef utilMesh2

% utilMesh2
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
        function obj = utilMesh2(val)
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
            %
            Mesh.Cn = [Mesh1.Cn; Mesh2.Cn+N];
            Mesh.CoorN = [Mesh1.CoorN; Mesh2.CoorN];
            Mesh.CoorV = [Mesh1.CoorV; Mesh2.CoorV];
            Mesh.Nsd = [Mesh1.Nsd; Mesh2.Nsd];
        end
        %
        % Extract a part of mesh
        function Mesh1 = ExtractMesh(obj,Pe)
            %
            Mesh = [obj.mesh];
            Mesh1 = Mesh;
            %
            Mesh1.Nsd = Mesh.Nsd(Pe);
            %if size(Mesh.CoorN,2) == 2, Mesh1.Surf = Mesh.Surf(Pe); end
            %
            Pn = unique(Mesh.Cn(Pe,:));
            %
            Cn = Mesh.Cn(Pe,:);
            %
            Mesh1.CoorN =  Mesh.CoorN(Pn,:);
            P1 = zeros(1,size(Mesh.CoorN,1));
            for k = 1:length(Pn), P1(Pn(k)) = k; end  % compresser la numérotation des noeuds
            Mesh1.Cn = P1(Cn);
            
            %
            
            if size(Mesh.CoorN,2) == 3
                Mesh1.CoorV=  Mesh.CoorV(Pe,:);
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
                    Mesh1 = ExtractMesh(utilMesh2(Mesh0(k)),Pe);
                    Mesh1 = MoveMesh(utilMesh2(Mesh1),[-xmax -ymax 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)>0 & Mesh0(k).CoorV(:,2)<0);
                    Mesh2 = ExtractMesh(utilMesh2(Mesh0(k)),Pe);
                    Mesh2 = MoveMesh(utilMesh2(Mesh2),[-xmax -ymin 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)<0 & Mesh0(k).CoorV(:,2)<0);
                    Mesh3 = ExtractMesh(utilMesh2(Mesh0(k)),Pe);
                    Mesh3 = MoveMesh(utilMesh2(Mesh3),[-xmin -ymin 0]);
                    %
                    Pe = find(Mesh0(k).CoorV(:,1)<0 & Mesh0(k).CoorV(:,2)>0);
                    Mesh4 = ExtractMesh(utilMesh2(Mesh0(k)),Pe);
                    Mesh4 = MoveMesh(utilMesh2(Mesh4),[-xmin -ymax 0]);
                    %
                    if length(Mesh0) == 1
                        Mesh = AssembMesh(utilMesh2(AssembMesh(utilMesh2(Mesh1),utilMesh2(Mesh2))),...
                                   utilMesh2(AssembMesh(utilMesh2(Mesh3),utilMesh2(Mesh4))));
                    else
                        Mesh(k) = AssembMesh(utilMesh2(AssembMesh(utilMesh2(Mesh1),utilMesh2(Mesh2))),...
                                   utilMesh2(AssembMesh(utilMesh2(Mesh3),utilMesh2(Mesh4))));
                    end
                end
            end
        
        end

    end
end

