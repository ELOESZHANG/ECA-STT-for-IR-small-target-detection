function [Pdist,C,S]=NL_wdist(f,patch_size,win_size,varargin)
% function NL_wdist finds the weighted distance in a gray/color image
% between patches of size patch_size and in a search window of size
% win_size. It returns only the K smallest distances between patches for 
% each pixel coordinate. The center of these patches are given in C.

% isgrad: If set to true then patch distances are computed based on the
% gradient of the input image.
% bc: type of extension for the image boundaries.
% ['symmetric'|'circular'|'zero'].

[bc,isgrad,weights,K]=process_options(varargin,'bc','symmetric',...
  'isgrad',false,'weights',ones(patch_size),'K',prod(win_size));

if K > prod(win_size)
  error('NL_wdist: The number of stored weights must be less than %d.',prod(win_size));
end

if any(mod(patch_size,2)-1)
  error('NL_wdist: The dimensions of the patch should be odd-sized.');
end

if any(mod(win_size,2)-1)
  error('NL_wdist: The dimensions of the search window should be odd-sized.');
end

if isscalar(weights)
  G=fspecial('gaussian',patch_size,weights);
else
  G=weights;
end

K=K-1;
[Nx, Ny, ~]=size(f);

Pdist=zeros([Nx Ny K]);
Cx=zeros([Nx Ny K]); % coordinates of the patch with the given similarity weights
Cy=zeros([Nx Ny K]);


T=reshape(1:Nx*Ny,[size(f,1),size(f,2)]);
% matrix which keeps the indexing of the image pixels.

[Y,X]=meshgrid(1:Ny,1:Nx);

if isgrad
  f=GradOp2D(f,bc);
end

ctr=1;
for kx=-floor(win_size(1)/2):floor(win_size(1)/2)
  for ky=-floor(win_size(2)/2):floor(win_size(2)/2)
    
    if ~(kx==0 && ky==0) % We do not check the (0,0)-shift case since in
      %this case each patch is compared to itself and the distance would be
      %zero. We add this weight at the end. (This is why we redefine above
      %K=K-1.)
      
      
      fs=shift(f,-[kx ky],bc);% fs(m,n)=f(m+kx,n+ky)
      e=abs(f-fs).^2;
      %e=(f-fs).^2;
      if isequal(bc,'zero')
        D=imfilter(e,G,'conv');
      else
        D=imfilter(e,G,'conv',bc);
      end
        
      if isgrad
        D=sum(sum(D,4),3);
      else
        D=sum(D,3);
      end
      %     fh=figure(ctr);
      %     msg=['kx=' num2str(kx,'%2.2f') ', ky=' num2str(ky,'%2.2f')];
      %     set(fh,'name',msg);
      %     imshow(D,[]); pause;
      
      if ctr <= K
        Pdist(:,:,ctr)=D;
        Cx(:,:,ctr)=repmat(kx,[Nx Ny 1 ]);
        Cy(:,:,ctr)=repmat(ky,[Nx Ny 1 ]);
      else
        % Check if the maximum element of Pdist is greater than the value of D
        % If yes, we keep the minimum weight. At the end we will have kept
        % the K smallest weights.
        [M,idx]=max(Pdist,[],3);
        idx=T+(idx-1)*Nx*Ny;
        % idx keeps the original location in Pdist of the elements stored in M,
        % i.e., M=Pdist(idx).
        mask=M > D;
        M=(1-mask).*M+mask.*D; % M=min(M,D)
        Pdist(idx)=M;
        idx(mask==0)=[]; % We store the new coordinates only for those
        % positions where M > D.
        Cx(idx)=kx;
        Cy(idx)=ky;
      end
      ctr=ctr+1;
      
    end
    
  end
end

% We add the weight for the (0,0)-shift case and the corresponding
% coordinates.

Pdist=cat(3,zeros(Nx,Ny),Pdist);
Cx=repmat(X,[1 1 K])+Cx;
Cx=cat(3,X,Cx);
Cy=repmat(Y,[1 1 K])+Cy;
Cy=cat(3,Y,Cy);

%We also need to fix Cx,Cy so as not to contain coordinates outside of the
%image domain Pdist=[1 Nx]x[1 Ny]

switch bc
  case 'symmetric'
    % pixels with negative or zero coordinates must be mapped to image
    % coordinates that belong to the image domain. For example if we
    % consider symmetric boundaries then a pixel (-m,-n) corresponds to the
    % pixel (m+1,n+1).
    Cx(Cx<=0)=abs(Cx(Cx<=0))+1; % Left boundaries that need to be fixed.
    Cx(Cx>=Nx+1)= 2*Nx+1-Cx(Cx>=Nx+1); % Right boundaries that need to be fixed.
    Cy(Cy<=0)=abs(Cy(Cy<=0))+1; % Upper boundaries that need to be fixed.
    Cy(Cy>=Ny+1)= 2*Ny+1-Cy(Cy>=Ny+1); % Lower boundaries that need to be fixed.
    
  case 'circular'
    % pixels with negative or zero coordinates must be mapped to image
    % coordinates that belong to the image domain. For example if we
    % consider circular boundaries then a pixel (-m,-n) corresponds to the
    % pixel (-m+N,n+1).
    Cx(Cx<=0)=Cx(Cx<=0)+Nx;
    Cx(Cx>=Nx+1)=Cx(Cx>=Nx+1)-Nx;
    Cy(Cy<=0)=Cy(Cy<=0)+Ny;
    Cy(Cy>=Ny+1)=Cy(Cy>=Ny+1)-Ny;
end

C=Cx+(Cy-1)*Nx; % scalar coordinate index.

if nargout > 2
  S=sort(sqrt(Pdist),3);
  S=S(:,:,floor(size(Pdist,3)/2));% distance of pixel (nx,ny) from neighbor (nx+k1,ny+k2)
  S(S==0)=eps;
  % S=sort(sqrt(Pdist),3,'ascend');
  % S=S(:,:,floor(K/2));% distance of pixel (nx,ny) from neighbor (nx+k1,ny+k2)
  % which will be used as the standard deviation of the produced weights.
  S=S(C);
  S=S.*repmat(S(:,:,1),[1 1 size(S,3)]);
end

function Df=GradOp2D(f,bc) %Gradient operator with boundary conditions bc

[r,c,nc]=size(f);
Df=zeros(r,c,2,nc);
Df(:,:,1,:)=shift(f,[-1,0],bc)-f; %f(i+1,j)-f(i,j)
Df(:,:,2,:)=shift(f,[0,-1],bc)-f; %f(i,j+1)-f(i,j)
