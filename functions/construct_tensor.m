function [D1,kk] = construct_tensor(D,numseq,i)
exlen = 2;
if numseq == 1
   D1 = D;
else
     if i ==1
        for ii = 1:(exlen+1)
            D1(:,:,ii) = D(:,:,i+ii-1);
        end
     elseif i == numseq
         for ii = 1:(exlen+1)
             D1(:,:,ii) = D(:,:,i+ii-(exlen+1));
         end
     else
         for ii = 1:(2*exlen+1)
             D1(:,:,ii) = D(:,:,i+ii-(exlen+1));
         end
     end
end
[~,~,kk]=size(D1);
