% 2020-11-06
% This matlab code implements the ECA-STT model for infrared small target
% detection
% Author: Lingyi Zhang

clc;
clear;
close all;

addpath('functions/')
datapath = 'images/seq1/';
imgDir = dir([datapath '*.bmp']);
saveDir = 'results/seq1/';

%% initial parameters
 lambda1 = 0.009;  % tuning
 L = 5.0;  % L, tuning 
 beta = 0.1;       %    tuning 
 t = 3;            % tuning
 
 lambda3 = 100;
 mu = 0.001;
 len = length(imgDir);
%% input IR image/sequence
 for i = 1:len
     img=imread(strcat(datapath,[num2str(i),'.bmp'])); 
     nn = ndims(img);
      if nn==3
          img= rgb2gray(img);
      end
     [m,n] = size(img);
     D(:,:,i) = double(img);
 end
%% low-rank and sparse recovery
for i = 1:t:len
     %% construct spatial-temporal tensor
     kk = min(t, len-i+1);
     D1 = D(:,:,i:(i+kk-1));
     lambda2 = L / sqrt(min(m,n)*kk);
     Wetensor = zeros(size(D1));
     %% compute target enhancement weight
      for ii = 1:kk
        img = D1(:,:,ii);
        [lambda11, lambda22] = structure_tensor_lambda(img,3);
        reCI = (lambda11+lambda22)./(lambda11.*lambda22);          % reciprocal of coner awareness indicator
        reEI = 1./(lambda11-lambda22);                             % reciprocal of edge awareness indicator
        WeTE = ((1+beta.^2).*reCI.*reEI)./((beta.^2.*reCI)+reEI);  % target enhancement weight
        Wetensor(:,:,ii)=WeTE;
      end
     %% compute tenNLTV weight 
        % Reference:
        % [1] W. Cao, J. Yao, J. Sun, and G. Han, ¡°A tensor-based nonlocal total variation model 
        % for multi-channel image recovery,¡± Signal Process., vol. 153, pp. 321¨C335, Dec. 2018.
      winsize = [11 11]; patchsize = [7 7];
      [D_STV_NL,C_STV_NL]=NL_wdist(D1,patchsize,winsize,'K',9,'bc','symmetric','isgrad',false);
      W_STV_NL=exp(-D_STV_NL/(0.25)^2);
     %% optimization via admm
      % fast solver: using PSVT operator to solve the tensor truncated nuclear norm
      % Reference:
      % [2] T. Oh, H. Kim, Y. Tai, J. Bazin, and I. Kweon, ¡°Partial sum minimization of singular values in rpca 
      % for low-level vision,¡± in 2013 IEEE Int. Conf.Comput. Vis. (ICCV), Los Alamitos, CA, USA, dec 2013, pp. 145¨C152.
%       [tenB, tenT ,tenN] = ECASTT_PSVT(D1,lambda1,lambda2,lambda3,mu,Wetensor, W_STV_NL,C_STV_NL);
      
      % original solver: using two-step iterative method to solve the tensor truncated nuclear norm
      % Reference:
      % S. Xue, W. Qiu, F. Liu, and X. Jin,¡°Truncated nuclear normregularization for low-rank tensor completion,¡± 2019. 
      % [Online].Available: https://arxiv.org/abs/1901.01997
      [tenB, tenT ,tenN] = ECASTT_TwoStep(D1,lambda1,lambda2,lambda3,mu,Wetensor, W_STV_NL,C_STV_NL);
     %% output the separated target image of the current frame
      for ii = 1:kk
         T_out = tenT(:,:,ii);
         T_out=im2uint8(mat2gray(T_out));
         figure; imshow(T_out);
         imwrite(T_out, [saveDir, num2str(i+ii-1),'.bmp']);
      end
end