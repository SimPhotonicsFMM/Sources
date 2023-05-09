function MatS = ProdMatS(varargin)

% ProdMatS 
%   S-matrix product 
%
% Syntax
%   MatS = ProdMatS(MatS1);  % MatS1 : Product of Array of S-matrices 
%   MatS = ProdMatS(MatS1,N);  % Product of MatS1 N times
%   MatS = ProdMatS(MatS1,MatS2,N);  % Product of (MatS1*MatS2) N times
%   MatS = ProdMatS(MatS1,MatS2,MatS3,...); %Product of a list of S-Matrix

% Date of the latest version : 13 February 2023
% Author : Mondher Besbes (LCF / CNRS / IOGS)


MatS1 = varargin{1};
%
if size(varargin,2) == 1 
    N = 1;
    %   MatS = ProdMatS(MatS1);  % MatS1 : tableau de matrices S
else
    if isscalar(varargin{1,2})
        N = varargin{1,2};
        %   MatS = ProdMatS(MatS1,N);  % Produit du tableau MatS1 N fois
    else
        MatS2 = varargin{1,2};
        MatS = Prod2MatS(MatS1,MatS2);
        if size(varargin,2) == 3
            if isscalar(varargin{1,3})
                %   MatS = ProdMatS(MatS1,MatS2,N);
                MatS = Prod1MatS(MatS,varargin{1,3});
            else
                %   MatS = ProdMatS(MatS1,MatS2,MatS3,...);
                for k = 3:size(varargin,2) 
                    MatS = Prod2MatS(MatS,varargin{1,k});
                end
            end
        end
        return
    end
end

%   MatS = ProdMatS(MatS1,N);  % Produit du tableau MatS1 N fois
MatS = MatS1(1,:);
%
for kc = 2:size(MatS1,1), MatS = Prod2MatS(MatS,MatS1(kc,:)); end
%
if N > 1, MatS = Prod1MatS(MatS,N); end
%
end

%%
function MatS2 = Prod1MatS(MatS1,N)
%
% Produit optimisé de matrice S de N fois
%

%
MatS2 = [];
if N == 1, MatS2 = MatS1; return, end
%
while N >= 1
    if mod(N,2)==1
        if isempty(MatS2)
            MatS2 = MatS1;
        else
            MatS2 = Prod2MatS(MatS1,MatS2);
        end
    end
    N = floor(N/2);
    if N > 0, MatS1 = Prod2MatS(MatS1,MatS1); end
end

end

%%
function MatS3 = Prod2MatS(MatS1,MatS2)
%
% Produit de deux matrices S
%
[n1,m1] = size(MatS1{1});
[n2,m2] = size(MatS2{1});

if length(MatS1) >= 3 
    N = size(MatS1{6},1); 
elseif length(MatS2) >= 3 
    N = size(MatS2{6},1);
else
    N = max([n1,m1,n2,m2])/2;
end
%
[n11,n12,n13,n14] = deal(m1-N,N,N,n1-N);
[n21,n22,n23,n24] = deal(N,m2-N,n2-N,N);
%   
if size(MatS1,2) >= 8
    MatS3 = MatS1;
elseif size(MatS2,2) >= 8
    MatS3 = MatS2;
end
%
VectS1 = full(MatS1{2});
VectS2 = full(MatS2{2});
MatS1 = MatS1{1};
MatS2 = MatS2{1};

%
if nnz(MatS1)/size(MatS1,1) < 10 && nnz(MatS2)/size(MatS2,1) <10
    
%MatS1 = sparse(MatS1);
%MatS2 = sparse(MatS2);
%

%tic
MatS3{1,1} = [sparse(n23,n11), MatS2(1:n23,n21+1:end); MatS1(n13+1:end,1:n11), sparse(n14,n22)] + ...
        ([MatS2(1:n23,1:n21) , sparse(n23,n12); sparse(n14,n21), MatS1(n13+1:end,n11+1:end)]*...
         ([speye(n13,n21) , -MatS1(1:n13,n11+1:end); -MatS2(n23+1:end,1:n21), speye(n24,n12)]\...
         [MatS1(1:n13,1:n11) , sparse(n13,n22); sparse(n24,n11), MatS2(n23+1:end,n21+1:end)]));
%toc
%MatS3{1,1} = full(MatS3{1,1});

else
% MatS1 = full(MatS1);
% MatS2 = full(MatS2);
%
%tic
% MatS3{1,1} = [zeros(n23,n11), MatS2(1:n23,n21+1:end); MatS1(n13+1:end,1:n11), zeros(n14,n22)] + ...
%         ([MatS2(1:n23,1:n21) , zeros(n23,n12); zeros(n14,n21), MatS1(n13+1:end,n11+1:end)]*...
%          ([eye(n13,n21) , -MatS1(1:n13,n11+1:end); -MatS2(n23+1:end,1:n21), eye(n24,n12)]\...
%          [MatS1(1:n13,1:n11) , zeros(n13,n22); zeros(n24,n11), MatS2(n23+1:end,n21+1:end)]));
% 
%toc
S2_11 = MatS2(1:n23,1:n21);
S2_12 = MatS2(1:n23,n21+1:end);
S2_21 = MatS2(n23+1:end,1:n21);
S2_22 = MatS2(n23+1:end,n21+1:end);
S1_11 = MatS1(1:n13,1:n11);
S1_12 = MatS1(1:n13,n11+1:end);
S1_21 = MatS1(n13+1:end,1:n11);
S1_22 = MatS1(n13+1:end,n11+1:end);

I = speye(n13,n21);
IS1 = I-S1_12*S2_21;
try dIS1 = decomposition(IS1,'lu'); catch, dIS1 = IS1; end
%InvIS1 = inv(IS1);

InvS1 = dIS1\S1_11;%InvIS1*S1_11;%
InvS2 = (dIS1\S1_12)*S2_22;%(InvIS1*S1_12)*S2_22;%
S2 = S1_22*S2_21;

MatS3{1,1} = [S2_11*InvS1 , S2_11*InvS2+S2_12 ; S1_21+S2*InvS1 , S1_22*S2_22+S2*InvS2];

%norm(MatS3{1,1}(:)-S(:))


end
%
if all(VectS1) == 0 && all(VectS2) == 0, MatS3{1,2} = VectS1; return, end
%
MatS3{1,2} = [VectS2(1:n23); VectS1(n13+1:end)] + ...
           [MatS2(1:n23,1:n21) , zeros(n23,n12); zeros(n14,n21), MatS1(n13+1:end,n11+1:end)]*...
           ([eye(n13,n21) , -MatS1(1:n13,n11+1:end); -MatS2(n23+1:end,1:n21), eye(n24,n12)]\...
           [VectS1(1:n13); VectS2(n23+1:end)]);

end
 





