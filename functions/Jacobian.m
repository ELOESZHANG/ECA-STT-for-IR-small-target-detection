function D = Jacobian(u,W,C,bc)
% Operator of the discrete Non Local Patch-based Jacobian
%
% u: Nx x Ny x Nc array of a vector-valued image with Nc-channels, defined
% on a Nx x Ny pixel grid
%
% W: Nx x Ny x Nw array with the Nw weights of each patch centered at
% location l=(nx,ny), 1<= nx <= Nx, 1<= ny <= Ny. These weights are inverse
% proportional to the patch distance between the two patches and are
% computed by the function NL_weights.
%
% C: Nx x Ny x Nw array with the center coordinates of the Nw patches
% which are related to the weights. 
% 
% bc: boundary condition type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
% D: 5D array with dimensions (Nx,Ny,2,Nw,Nc), which for each pixel (i,j)
% contains the 2 x Nw x Nc Non Local Patch-based Jacobian D(i,j,:,:,:). 
if nargin < 4
  bc='symmetric';
end

[Nx,Ny,Nc] = size(u);
% P = Nx*Ny;

Nw=size(W,3);

% if ~(mod(Nw, 2)) % if not all [NGx,NGy] are odd numbers
%     error('The number of weights must be odd.');
% end

grad_c = zeros(Nx,Ny,2,Nc);%2->d,[nX,nY]->size of image,Nc->channel numbers
for i_chan=1:Nc
    % gradient of each image channel:
    grad_c(:,:,:,i_chan) = GradOp2D(u(:,:,i_chan),bc);
end

%grad_c = VGradOp2D(u,bc);
clear u;


D = zeros(Nx,Ny,2,Nw,Nc);%2*(Nw*Nc)->d*(B*L)
D(:,:,:,1,1:Nc) = grad_c.*repmat(W(:,:,1),[1 1 2 Nc]);
idx=reshape(kron((0:2*Nc-1)*Nx*Ny,ones(Nx,Ny)),[Nx Ny 2 Nc]);
for k=2:Nw
    idx2=repmat(C(:,:,k), [1 1 2 Nc])+idx;
    T=grad_c(idx2).*repmat(W(:,:,k),[ 1 1 2 Nc]);   
    D(:,:,:,k,1:Nc) = T;
end

function Df=GradOp2D(f,bc) %Gradient operator with boundary conditions bc

[r,c]=size(f);
Df=zeros(r,c,2);
Df(:,:,1)=shift(f,[-1,0],bc)-f; %f(i+1,j)-f(i,j)
Df(:,:,2)=shift(f,[0,-1],bc)-f; %f(i,j+1)-f(i,j)


% function Df=VGradOp2D(f,bc) %Gradient operator with boundary conditions bc
% 
% [r,c,k]=size(f);
% Df=zeros(r,c,2,k);
% Df(:,:,1,:)=shift(f,[-1,0,0],bc)-f; %f(i+1,j,k)-f(i,j,k)
% Df(:,:,2,:)=shift(f,[0,-1,0],bc)-f; %f(i,j+1,k)-f(i,j,k)
