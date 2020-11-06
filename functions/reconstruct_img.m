function [E_out] = reconstruct_img(tenT,numseq,i)
exlen = 2;
if numseq ==1
    E_out = tenT;
else
     if i == 1
        E_out=tenT(:,:,1);
     elseif i == numseq
        E_out=tenT(:,:,exlen+1);
     else
        E_out=tenT(:,:,exlen+1);
     end 
end