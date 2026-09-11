function [x0,y0,z0] = sf3d(n0,a0,angle0, Dep0, nu,nv)

%
u = linspace(-pi,pi,nu);
v = linspace(-pi/2,pi/2,nv);
nu = length(u);
nv = length(v);
%
if mod(nv,2) == 0, warning('Choose odd number for angle Phi'); end
%
%[x,y,z] = deal(nan(nu,nv));
[x0,y0,z0] = deal(nan(size(n0,1),nu,nv));
%
a1 = [1 1];
%
for k = 1:size(n0,1)
    [n,a,angle,Dep] = deal(n0(k,:),a0(k,:),angle0(k),Dep0(k,:));

    %
    for i = 1:nu
    for j = 1:nv
      raux1 = abs(1 / a1(1) * abs(cos(n(1) .* u(i) / 4))) .^ n(3) + abs(1 / a1(2) * abs(sin(n(1) * u(i) / 4))) .^ n(4);
      r1 = abs(raux1) .^ (- 1 / n(2));
      raux2 = abs(1 / a1(1) * abs(cos(n(1) * v(j) / 4))) .^ n(3) + abs(1 / a1(2) * abs(sin(n(1) * v(j) / 4))) .^ n(4);
      r2 = abs(raux2) .^ (- 1 / n(2));
      x(i, j) = r1 * cos(u(i)) * r2 * cos(v(j));
      y(i, j) = r1 * sin(u(i)) * r2 * cos(v(j));
      z(i, j) = r2 * sin(v(j));
    end
    end
    %
    x = x-mean([max(x(:)), min(x(:))]);
    y = y-mean([max(y(:)), min(y(:))]);
    z = z-mean([max(z(:)), min(z(:))]);
    scale = [2*a(1)/(max(x(:))-min(x(:))) 2*a(2)/(max(y(:))-min(y(:))) 2*a(3)/(max(z(:))-min(z(:)))];
    [x,y,z] = deal(x*scale(1) , y*scale(2) , z*scale(3));
    xyz = [1 0 0; 0 cos(angle) sin(angle);0 -sin(angle) cos(angle)]*[x(:)';y(:)';z(:)'];
    x = reshape(xyz(1,:),size(x));
    y = reshape(xyz(2,:),size(y));
    z = reshape(xyz(3,:),size(z));
    %
    [x0(k,:,:),y0(k,:,:),z0(k,:,:)] = deal(x , y , z+Dep(3) );
    %
end
%
%if size(n0,1) == 1, [x0,y0,z0] = deal(squeeze(x0),squeeze(y0),squeeze(z0)); end
%
figure, hold on
% if size(n0,1) == 1
%     mesh(x0,y0,z0); axis equal    
% else
    for k = 1:size(n0,1)
        Dep = Dep0(k,:);
        mesh(squeeze(x0(k,:,:) + Dep(1)), squeeze(y0(k,:,:) + Dep(2)), squeeze(z0(k,:,:))); axis equal, axis tight
    end
%end

end
