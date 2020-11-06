function v = AdjJacobian(D,W,idx,n,I,bc)
% Computes the adjoint of the discrete Non Local Patch-based Jacobian
% (see JacobianOp2D_NL).
% D: 5D array with dimensions (Nx,Ny,2,Nw,Nc), which for each pixel (i,j)
% contains the 2 x Nw x Nc Non Local Patch-based Jacobian D(i,j,:,:,:). 
% W: Nx x Ny x Nw array with the Nw weights of each patch centered at
% location l=(nx,ny), 1<= nx <= Nx, 1<= ny <= Ny. These weights are inverse
% proportional to the patch distance between the two patches and are
% computed by the function NL_weights.
%
% C: Nx x Ny x Nw array with the center coordinates of the Nw patches
% which are related to the weights. From C we obtain the following three
% input parameters using function compute_Nsum_indices():
% idx: vector of size (Nx x Ny x Nw) x 1 ( [~,idx]=sort(C(:)); )
% n  : vector of size (Nx x Ny) x 1 ( n=histc(C(:),1:Nx*Ny); )
% I  : vector of size (Nx x Ny) x 1 ( I=cumsum(n)'; I=[1 I(1:end-1)+1];)
%
% bc: boundary conditions type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
%
% v: Nx x Ny x Nc array (same dims as vectorial image) with the result of
% the adjoint operator on D
if nargin < 6
  bc='symmetric';
end


[Nx,Ny,Nw] = size(W);
Nc=size(D,5);

F = zeros(Nx,Ny,2,Nw,Nc);

v=zeros(Nx,Ny,Nc);
tmp=zeros(Nx,Ny,2);
for k=1:Nc
  for m=1:Nw
    F(:,:,:,m,k)=(repmat(W(:,:,m),[1 1 2])).*D(:,:,:,m,k);
  end
  
  tmp(:,:,1)=Nsum_mex(F(:,:,1,1:Nw,k),idx-1,n,I-1);
  tmp(:,:,2)=Nsum_mex(F(:,:,2,1:Nw,k),idx-1,n,I-1);
  % We need to substract one from idx and I since Nsum_mex is a c file and 
  % therefore the array indexing starts from 0 instead of 1. 
  v(:,:,k)=AdjGradOp2D(tmp,bc);  
end


function g=AdjGradOp2D(P,bc) %Adjoint gradient operator (i.e. -div)

P1=P(:,:,1);
P1=shiftAdjST(P1,[-1,0],bc)-P1; % P1(i-1,j)-P1(i,j)
P2=P(:,:,2);
P2=shiftAdjST(P2,[0,-1],bc)-P2; % P2(i,j-1)-P2(i,j)
g=P1+P2;