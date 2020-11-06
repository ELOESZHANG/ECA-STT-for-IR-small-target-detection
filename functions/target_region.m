function [T1,label] = target_region(T)
% [m,n] = size(T);
[r,c] = find(T~=0);
len = length(r);
label = zeros(len,1);
for i = 1:len
    k = max(label(:));
    label(i)=k+1;
    for j = 1:i
        if abs(r(j)-r(i))<=1 && abs(c(j)-c(i))<=1
           if label(i) >label(j)
              label(i)=label(j);
           else
               label(j)=label(i);
           end

        end
    end
end
areanum = max(label(:));
areas =zeros(areanum,1);
pixelnum=zeros(areanum,1);
for ii = 1:areanum
    for i=1:len
        if label(i)==ii
            areas(ii)=areas(ii)+T(r(i),c(i))/255;
            pixelnum(ii)=pixelnum(ii)+1;
        end
    end
    if pixelnum(ii)~=0
        areas(ii)=areas(ii)./pixelnum(ii);
    end
end

% areas=areas./pixelnum
[x,y]=find(areas==max(areas(:)));
targetlabel = x;
T1 = T;
for i=1:len
    if label(i)~=targetlabel
      T1(r(i),c(i))=0;
    end
end

        