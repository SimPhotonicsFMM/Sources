function [x0,y0] = sf2d(varargin)

% sf2d 
%   Generate geometric shape from J. Gielis superformula
%   https://en.wikipedia.org/wiki/Superformula
%
% Syntax
%   [x0,y0] = sf2d(n0, a0, angle0, Dep0, Np);               
%   [x0,y0] = sf2d(n0, a0, angle0, Dep0);               
%   [x0,y0] = sf2d(n0, a0, angle0);               
%   [x0,y0] = sf2d(n0, a0);               
%   [x0,y0] = sf2d(Geom);    % see SetGeom           
%
% Description
%   n0 : List of parameters [m, n1, n2, n3]
%   a0 : Radius [a b]
%   angle0 : Angle of rotation (rd)
%   Dep0   : Translation in x and y
%   Np     : Nomber of points (Np = 401 by default)
%   Geom   : Structure with different parameters (see SetGeom)
%            [ mn , ab et/ou (lix, liy) , Angle , Dep , Np ]
%
%   x0 : Abscissa of points
%   y0 : Ordonates of points
%
% Example
%   mn = [m, n1, n2, n3]  ; ab = [rx ry]
%   [4 1 1 1] : Diamond
%   [4 2 2 2] : Circle or ellipse
%   [4 n n n] : n>2 rounded square or rectangle 
%   [3 n 2*n 2*n] : rounded triangle n>1
%   mn = 4 : Rectangle or square 
%
%   [xv,yv] = sf2d([3 4 8 8 ; 3 4 8 8], [1 1; 1 1], ...
%                [pi 0], [-1.1 0 ; 1.1 0]);
%   figure, plot(xv',yv'), axis equal

% Date of the latest version : 10 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)


if isstruct(varargin{1,1})
    %
    Geom = varargin{1,1};
    n0 = Geom.mn;
    if isfield(Geom,'ab'), a0 = Geom.ab; else, a0 = [Geom.lix Geom.liy]/2; end
    if isfield(Geom,'Np'), Np = Geom.Np; else, Np = 181; end
    if isfield(Geom,'Dep'), Dep0 = Geom.Dep; else, Dep0 = zeros(size(a0)); end
    if isfield(Geom,'Angle'), angle0 = Geom.Angle; else, angle0 = zeros(size(n0,1),1); end
    
    [x0,y0] = sf2d(n0, a0, angle0, Dep0, Np);

    return
end

[n0, a0 ] = deal(varargin{1,1},varargin{1,2});

if nargin < 5, Np = 401; else, Np = varargin{1,5}; end
if nargin < 4, Dep0 = zeros(size(a0)); else, Dep0 = varargin{1,4}; end
if nargin < 3, angle0 = zeros(size(n0,1),1); else, angle0 = varargin{1,3}; end


u = linspace(0,2*pi,Np);
if size(n0,2) == 4
    [x0,y0] = deal(nan(size(n0,1),length(u)));
else
    [x0,y0] = deal(nan(size(n0,1),5));  % quadrangle
end

%
for k = 1:size(n0,1)
    %
    [n,a,angle,Dep] = deal(n0(k,:),a0(k,:),angle0(k),Dep0(k,:));
    if length(n) == 4
        a1 = [1 1];
        raux = abs(1 / a1(1) .* abs(cos(n(1) * u / 4))) .^ n(3) + ...
               abs(1 / a1(2) .* abs(sin(n(1) * u / 4))) .^ n(4);
        %
        r = abs(raux) .^ (- 1 / n(2));
        x = r .* cos(u);
        y = r .* sin(u);
    else
        x = [-a(1)  a(1) a(1) -a(1) -a(1)];
        y = [-a(2)  -a(2) a(2) a(2) -a(2)];
    end
    x = x-mean([max(x), min(x)]);
    y = y-mean([max(y), min(y)]);
    scale = [2*a(1)/(max(x)-min(x)) 2*a(2)/(max(y)-min(y))];
    [x,y] = deal(x*scale(1) , y*scale(2) );

    xy = [cos(angle) -sin(angle); +sin(angle) cos(angle)]*[x;y];
    [x0(k,:),y0(k,:)] = deal(xy(1,:) + Dep(1) , xy(2,:) + Dep(2));
    %
end

end