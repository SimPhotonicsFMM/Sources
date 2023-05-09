function VisuContour(Mesh,NumPlan,LineColor)

% VisuContour 
%   Plot contours of subdomains of a 2D or 3D mesh
%
% Syntax
%   VisuContour(Mesh,LineColor)
%   VisuContour(Mesh)
%
% Description
%   Mesh : Mesh of a structure
%   LineColor: Color of lines 
%

% Date of the latest version : 14 March 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)

if nargin < 2, NumPlan = 0;  end
if nargin < 3, LineColor = 'k'; end

if length(Mesh)>1
    for k = 1:length(Mesh)
        VisuContour(Mesh(k),NumPlan,LineColor);
    end
    return
end


Flag = zeros(length(Mesh.ExtAr),1);

for ia = 1:length(Mesh.ExtAr)
    [Pe,~] = find(ismember(abs(Mesh.Ca),ia));
    if numel(Pe) == 1 
        Flag(ia) = 1; 
    elseif numel(Pe) >= 2
        if Mesh.Nsd(Pe(1)) ~= Mesh.Nsd(Pe(2))
            Flag(ia) = 2;
        end
    end
end

% -------------------------- Sous-domaines -------------------------%
x = Mesh.CoorN(:,1);
y = Mesh.CoorN(:,2);
if size(Mesh.CoorN,2) == 3, z = Mesh.CoorN(:,3); end

%
hold on
Pa = find(Flag == 1 | Flag == 2); 
for k = 1:1:length(Pa) 
    ia = Pa(k);
    switch NumPlan
        case 0
            if size(Mesh.CoorN,2) == 3
                line([x(Mesh.ExtAr(ia,1)) x(Mesh.ExtAr(ia,2))],...
                [y(Mesh.ExtAr(ia,1)) y(Mesh.ExtAr(ia,2))],...
                [z(Mesh.ExtAr(ia,1)) z(Mesh.ExtAr(ia,2))],...
                'Color',LineColor,'LineWidth',1), 
            else
                line([x(Mesh.ExtAr(ia,1)) x(Mesh.ExtAr(ia,2))],...
                [y(Mesh.ExtAr(ia,1)) y(Mesh.ExtAr(ia,2))],...
                'Color',LineColor,'LineWidth',1),
            end
        case 1
                line([y(Mesh.ExtAr(ia,1)) y(Mesh.ExtAr(ia,2))],...
                [z(Mesh.ExtAr(ia,1)) z(Mesh.ExtAr(ia,2))],...
                'Color',LineColor,'LineWidth',1),
        case 2
                line([x(Mesh.ExtAr(ia,1)) x(Mesh.ExtAr(ia,2))],...
                [z(Mesh.ExtAr(ia,1)) z(Mesh.ExtAr(ia,2))],...
                'Color',LineColor,'LineWidth',1),
        case 3
                line([x(Mesh.ExtAr(ia,1)) x(Mesh.ExtAr(ia,2))],...
                [y(Mesh.ExtAr(ia,1)) y(Mesh.ExtAr(ia,2))],...
                'Color',LineColor,'LineWidth',1),
    end
end

axis equal

return