function Mesh = MeshFreeForm(geom,MyFun,pas,alphaX,alphaY)

% pas : Step of Thickness of layer

% Initial mesh for homogenous mediun
geom0 = SetGeom(geom,'SaveOption',0);
Mesh0 = MeshLayer(geom0);
%
xn = unique(Mesh0.CoorN(:,1));%
yn = unique(Mesh0.CoorN(:,2));
[Xn,Yn] = ndgrid(xn,yn);
%
xp = unique(Mesh0.CoorV(:,1));
yp = unique(Mesh0.CoorV(:,2));
[X,Y] = ndgrid(xp,yp);
%
if nargin < 4,  alphaX = 1; alphaY = 1;end
if nargin < 5,  alphaY = 1;end
%
P = abs(X)<=alphaX*max(xp) & abs(Y)<=alphaY*max(xp);
%
% Add position min and max of each layer 
sz = struct();
[zmin,zmax] = deal(zeros(1,length(pas)));
%
figure, hold on
for k = 1:length(pas)
    sz(k).Z = zeros(size(X));
    sz(k).Z(P) = MyFun(X(P),Y(P));
    sz(k).Z = sz(k).Z +sum(pas(1:end-k+1));
    Zn = MyFun(Xn,Yn)+sum(pas(1:end-k+1));
    zmin(k) = min(Zn(:)); %min(sz(k).Z(:));
    zmax(k) = max(Zn(:)); %max(sz(k).Z(:));
    %zmin(k) = min(sz(k).Z(:));
    %zmax(k) = max(sz(k).Z(:));
    mesh(X,Y,sz(k).Z), axis equal; axis tight; view([-30 30]); drawnow
end
% for k = 1:length(pas)
%     sz(k).Z = MyFun(X,Y)+sum(pas(1:end-k+1));
%     Zn = MyFun(Xn,Yn)+sum(pas(1:end-k+1));
%     zmin(k) = min(Zn(:)); %min(sz(k).Z(:));
%     zmax(k) = max(Zn(:)); %max(sz(k).Z(:));
% end
%
% Update Mesh
geom = SetGeom(geom,'npz',[geom.npz(1) zmin zmax]);
Mesh0 = MeshLayer(geom);
Mesh0.Nsd(:) = max(geom.NumSD)+1;
%
% Update subdomain numbering
Mesh = Mesh0;
%
for k = 1:length(pas)
    Z = sz(k).Z;
    for kp = 1:length(X(:))
        % P = find(abs(Mesh0.CoorV(:,1)-X(kp))<1e-10 & ...
        %          abs(Mesh0.CoorV(:,2)-Y(kp))<1e-10);
        P = kp:length(xp)*length(yp):size(Mesh0.CoorV,1);
        Pz = find(Mesh.CoorV(P,3)<=Z(kp));
        Mesh.Nsd(P(Pz)) = geom.NumSD(k);
        % Pz = find(Mesh.CoorV(:,3)<=Z(kp));
        % Mesh.Nsd(Pz) = geom.NumSD(k);
    end
end
% Update facet numbering
Mesh = NumFacet(utilMesh(Mesh));
%
end