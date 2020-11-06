function [tenB, tenT,tenN] = ECASTT_TwoStep(tenD,lambda1,lambda2,lambda3,mu,wepara,W,C)
%% initialize
rho = 1.5;
tol = 10^(-7);
cg_tol=1e-4;cg_iter=10;
[h,w,t] = size(tenD);
d=2;L=9;B=t;   
R = min(h,w);  %%truncated rank
tenB = tenD;
tenT = zeros(size(tenD));
tenN = zeros(size(tenD));
tenY1 = zeros(size(tenD));tenY2 = zeros(size(tenD));tenY3 = zeros(h,w,d,L,B);tenY4 = zeros(size(tenD));
normD = norm(tenD(:));
reweightTen = ones(size(tenD)) .* wepara;
[idx,n,I]=compute_Nsum_indices(C);
bc='symmetric';

iter = 0;
converged = false;
while ~converged
    iter = iter + 1;
   %% Update Z
    if t >= 2
    tenZ = prox_tnn(tenB-tenY1/mu,1/mu);
    else
    [UMat, SMat, VMat] = svd(double((tenB-tenY1/mu)), 'econ');
    THSMat = max(SMat-1/mu, 0);
    tenZ=UMat*THSMat*VMat';
    end
   %% Update M
    tenB = tenZ;
    if t >= 2
    [U, ~, V] = t_SVD(tenB);
    else
    [U,~,V]=svd(tenB);
    end
    G = tran(U(:, 1:R, :)); H = tran(V(:, 1:R, :));
    GtH = tprod(tran(G), H);
    tenM = tenB+( GtH-tenY2) / mu;    
   %% Update L  
    JB = Jacobian(tenB,W,C,bc);
    tenL = softThres(JB-tenY3/mu,lambda1*ones(size(JB))/mu);
   %% Update B
    JL = AdjJacobian(tenL,W,idx,n,I,bc);
    JY3 = AdjJacobian(tenY3,W,idx,n,I,bc);
    temp1 = @(tenB)JTJ(tenB,W,idx,n,I,bc)+3*tenB;
    temp2 = (tenY1+tenY2+JY3+tenY4+mu*(tenZ+tenM+JL+tenD-tenT-tenN))/mu;
    tenB = cg_solver(temp1,temp2,tenB,cg_iter,cg_tol);
   %% Update T
    tenT = prox_non_neg_l1(tenD-tenB-tenN+tenY4/mu, lambda2*reweightTen/mu);
   %% Update W_ReTE
    reweightTen = 1 ./ (abs(tenT) + 0.01) .* wepara;
   %% Update N
    tenN = (mu*(tenD-tenB-tenT)+tenY4)/(mu+2*lambda3);
   %% Update Y
     tenY1 = tenY1+mu*(tenZ-tenB);
     tenY2 = tenY2+mu*(tenM-tenB);
     tenY3 = tenY3+mu*(tenL-Jacobian(tenB,W,C,bc));
     tenY4 = tenY4+mu*(tenD-tenB-tenT-tenN);
   %% Update mu
     mu = rho*mu;
   %% Output
      stopCriterion = norm(tenD(:) - tenB(:) - tenT(:) - tenN(:)) / normD; 
      if stopCriterion < tol
        converged = true;
     end
    disp(['#Iteration ' num2str(iter) ' |T|_0 ' ...
        num2str(sum(tenT(:) > 0)) ...
        ' stopCriterion ' num2str(stopCriterion)]);        
end