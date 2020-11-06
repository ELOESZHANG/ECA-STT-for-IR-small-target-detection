function v = JTJ(f,W,idx,n,I,bc)
% Computes the composite operation J^T(J(f)) where J is the discrete Non 
% Local Patch-based Jacobian(see JacobianOp2D_NL) and J^T is its adjoint 
% (see AdjJacobianOp2D_NL). Note that JTJOp2D_NL(f,sqrt(W),idx,n,I) gives
% exactly the same result with 
% AdjJacobianOp2D_NL(JacobianOp2D_NL(f,sqrt(W),C),sqrt(W),idx,n,I) but it 
% is much faster. 
% f: Nx x Ny x Nc array which represents a vector-valued image of Nc
% channels.
%
% W: Nx x Ny x Nw array with the Nw weights of each patch centered at
% location l=(nx,ny), 1<= nx <= Nx, 1<= ny <= Ny. These weights are inverse
% proportional to the patch distance between the two patches and are
% computed by the function NL_weights.
%
% C: Nx x Ny x Nw array with the center coordinates of the Nw patches
% which are related to the weights. From C we obtain the following three
% input parameters using function compute_Nsum_indices():
%
% idx: vector of size (Nx x Ny x Nw) x 1 ( [~,idx]=sort(C(:)); )
% n  : vector of size (Nx x Ny) x 1 ( n=histc(C(:),1:Nx*Ny); )
% I  : vector of size (Nx x Ny) x 1 ( I=cumsum(n)'; I=[1 I(1:end-1)+1];)
%
% bc: boundary conditions type: 'symmetric' |'circular'|'zero'.
%
% (OUTPUT):
%
% v: Nx x Ny x Nc array (same dims as vectorial image).
if nargin < 6
  bc='symmetric';
end

[Nx,Ny,Nc]=size(f);

W=Nsum_mex(W.^2,idx-1,n,I-1);

Wgrad = zeros(Nx,Ny,2,Nc);
v=zeros(Nx,Ny,Nc);
for i_chan=1:Nc
    % gradient of each image channel:
    Wgrad(:,:,:,i_chan) = repmat(W,[1 1 2]).*GradOp2D(f(:,:,i_chan),bc);
    v(:,:,i_chan)=AdjGradOp2D(Wgrad(:,:,:,i_chan),bc);
end


function g=AdjGradOp2D(P,bc) %Adjoint gradient operator (i.e. -div)

P1=P(:,:,1);
P1=shiftAdjST(P1,[-1,0],bc)-P1; % P1(i-1,j)-P1(i,j)
P2=P(:,:,2);
P2=shiftAdjST(P2,[0,-1],bc)-P2; % P2(i,j-1)-P2(i,j)
g=P1+P2;

function Df=GradOp2D(f,bc) %Gradient operator with boundary conditions bc

[r,c]=size(f);
Df=zeros(r,c,2);
Df(:,:,1)=shift(f,[-1,0],bc)-f; %f(i+1,j)-f(i,j)
Df(:,:,2)=shift(f,[0,-1],bc)-f; %f(i,j+1)-f(i,j)