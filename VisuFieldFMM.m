function VisuFieldFMM(Vect,X,Y,Z)

% VisuFieldFMM 
%   Plot the field "Vect" along an axis or on grid
%
% Syntax
%   VisuFieldFMM(Vect,X,Y,Z) % structure 3D
%   VisuFieldFMM(Vect,X,Y)   % structure 2D
%
% Description

%   Vect : Field vector
%
% Example
%   Run the example of the function "CulculFieldFMM"
%   figure, VisuFieldFMM(abs(VectE(:,3)),x,y,z)

% Date of the latest version : 31 January 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)


if nargin == 4
    if isvector(X) && isvector(Y) && isvector(Z)
        [x,y,z] = ndgrid(X,Y,Z);
    else
        [x,y,z] = deal(X,Y,Z);
    end
elseif nargin == 3
    if isvector(X) && isvector(Y)
        [x,y] = ndgrid(X,Y);
    else
        [x,y] = deal(X,Y);
    end
    z = y;
end
%
if isvector(squeeze(x)) && sum(diff(x(:))) ~= 0
    PlotCourbe(x(:),Vect,'x','','','r')
elseif isvector(squeeze(y)) && sum(diff(y(:))) ~= 0
    PlotCourbe(y(:),Vect,'y','','','r')
elseif isvector(squeeze(z)) && sum(diff(z(:))) ~= 0
    PlotCourbe(z(:),Vect,'z','','','r')
else
    
    if sum(diff(y(:))) == 0
        [N1,N2] = size(x);
        pcolor(squeeze(x),squeeze(z),reshape(Vect,N1,N2));
        title('XZ')
    elseif sum(diff(x(:))) == 0
        [~,N1,N2] = size(y);
        pcolor(squeeze(y),squeeze(z),reshape(Vect,N1,N2));
        title('YZ')
    elseif sum(diff(z(:))) == 0 || norm(y(:)-z(:))<eps
        [N1,N2] = size(x);
        pcolor(squeeze(x),squeeze(y),reshape(Vect,N1,N2));
        title('XY')
    end
    shading interp, axis equal, colormap jet; drawnow
end

end