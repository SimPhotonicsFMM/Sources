function [ht,mn,ab,DepXY,npx,npy] = GeomGrating2D(N,mn0,h,rb,rh,X0,Y0)
%
if nargin == 5, [X0,Y0] = deal(zeros(size(h))); end
if nargin == 6, Y0 = zeros(size(X0)); end
%
if numel(rb) == 1
    [rbx,rby] = deal(rb);
else
    rbx = rb(:,1)'; rby = rb(:,2)';
end
%
if numel(rh) == 1
    [rhx,rhy] = deal(rh);
else
    rhx = rh(:,1)'; rhy = rh(:,2)';
end
%
if all(X0(:)==0), nx0 = [30 2]; else, nx0 = ceil(size(h,1)*30); end
if all(Y0(:)==0), ny0 = [30 2]; else, ny0 = ceil(size(h,2)*30); end
%
Dx = max(X0(:))-min(X0(:));
DeltaX0 = Dx/nx0(1);
Dy = max(Y0(:))-min(Y0(:));
DeltaY0 = Dy/ny0(1);
%
ht = [];
ab = {};
DepXY = {};
mn = {};
npx = {};
npy = {};
%
h = h(:)';
N0 = N;
Nh = length(h);
X = X0(:); Y = Y0(:);
h10 = min(h)/N;
rhx0 = rhx;
rhy0 = rhy;
%
for kh = 1:Nh
    [h1,P] = min(h);
    rhx = rhx0+(h-h1).*(rbx-rhx0)./h;
    rhy = rhy0+(h-h1).*(rby-rhy0)./h;
    %
    rx = rhx'+(rbx-rhx)'*(1-(2*(1:N)-1)/(2*N));
    ry = rhy'+(rby-rhy)'*(1-(2*(1:N)-1)/(2*N));
    ht = [ h1/N*ones(1,N) ht]; % Hauteur des tranches
    %
    for k = 1:N 
        ab = [[rx(:,k) ry(:,k)] ab];
        if size(mn0,1) == length(h)
            mn = [mn0 mn]; 
        else
            mn = [repmat(mn0,length(h),1) mn]; 
        end
        %np = [[ceil(sqrt(length(h))*30) 2] np];
        %
        if mn0 == 4 
            nx = 2; ny = 2;
        else
            if DeltaX0 > min(rx(:,k))/2 
                DeltaX = min(rx(:,k))/4; nx = ceil(Dx/DeltaX);
            else
                DeltaX0 = Dx/nx0(1); nx = nx0;
            end
            %
            if DeltaY0 > min(ry(:,k))/2 
                DeltaY = min(ry(:,k))/4; ny = ceil(Dy/DeltaY);
            else
                DeltaY0 = Dy/ny0(1); ny = ny0;
            end

        end
        npx = [nx npx];
        npy = [ny npy];
        DepXY = [[X(:) Y(:)] DepXY];
    end 
    h = h-h1;
    h(P) = []; 
    rbx = rhx;
    rbx(P) = []; 
    if length(rhx0)>1, rhx0(P) = []; end
    rby = rhy;
    rby(P) = []; 
    if length(rhy0)>1, rhy0(P) = []; end
    N = ceil(min(h)/h10);
    %
    X(P) = []; Y(P) = []; 
    N0 = [N N0];
end

end