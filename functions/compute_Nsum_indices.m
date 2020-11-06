function [idx, n, I]=compute_Nsum_indices(C)

% idx: vector of size (Nx x Ny x Nw) x 1 ( [~,idx]=sort(C(:)); )
% n  : vector of size (Nx x Ny) x 1 ( n=histc(C(:),1:Nx*Ny); )
% I  : vector of size (Nx x Ny) x 1 ( I=cumsum(n)'; I=[1 I(1:end-1)+1];)

N=numel(C)/size(C,ndims(C));
[~,idx]=sort(C(:));
n=histc(C(:),1:N);
I=cumsum(n); I=[1; I(1:end-1)+1];
