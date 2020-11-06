function y=shiftAdjST(y,e,bc)


%dim_size=ndims(y);
ys=zeros(size(y));

% if dim_size ~= 4
%   error('shiftAdjST::Only 3D arrays can be processed');
% end

% if e > dim_size
%   e(dim_size+1:end)=[];
% elseif (e < dim_size)
%   e(dim_size)=0;
% end

switch bc
  case 'symmetric'
    if e(1) > 0
      ys(1:end-e(1),:,:,:) = y(1+e(1):end,:,:,:);
      ys(1:e(1),:,:,:) = ys(1:e(1),:,:,:) + y(e(1):-1:1,:,:,:);
    elseif e(1) < 0
      ys(1-e(1):end,:,:,:) = y(1:end+e(1),:,:,:);
      ys(end+e(1)+1:end,:,:,:) = ys(end+e(1)+1:end,:,:,:) + y(end:-1:end+e(1)+1,:,:,:);
    else
      ys=y;
    end
    
    y=zeros(size(ys));
    
    if e(2) > 0
      y(:,1:end-e(2),:,:) = ys(:,1+e(2):end,:,:);
      y(:,1:e(2),:,:) = y(:,1:e(2),:,:) + ys(:,e(2):-1:1,:,:);
    else
      y(:,1-e(2):end,:,:) = ys(:,1:end+e(2),:,:);
      y(:,end+e(2)+1:end,:,:) = y(:,end+e(2)+1:end,:,:) + ys(:,end:-1:end+e(2)+1,:,:);
    end
  case 'circular'
    if e(1) > 0
      ys(1:end-e(1),:,:,:)=y(1+e(1):end,:,:,:);
      ys(end-e(1)+1:end,:,:,:)=y(1:e(1),:,:,:);
    else
      ys(1-e(1):end,:,:,:)=y(1:end+e(1),:,:,:);
      ys(1:-e(1),:,:,:)=y(end+e(1)+1:end,:,:,:);
    end
    
    y=zeros(size(ys));
    
    if e(2) > 0
      y(:,1:end-e(2),:,:)=ys(:,1+e(2):end,:,:);
      y(:,end-e(2)+1:end,:,:)=ys(:,1:e(2),:,:);
    else
      y(:,1-e(2):end,:,:)=ys(:,1:end+e(2),:,:);
      y(:,1:-e(2),:,:)=ys(:,end+e(2)+1:end,:,:);
    end
  case 'zero'
    if e(1) > 0
      ys(1:end-e(1),:,:,:)=y(1+e(1):end,:,:,:);
    else
      ys(1-e(1):end,:,:,:)=y(1:end+e(1),:,:,:);
    end
    
    y=zeros(size(ys));
    
    if e(2) > 0
      y(:,1:end-e(2),:,:)=ys(:,1+e(2):end,:,:);
    else
      y(:,1-e(2):end,:,:)=ys(:,1:end+e(2),:,:);
    end 
  otherwise
      error('shift:InvalidShiftType','%s','Unknown boundary conditions.');    
    
end