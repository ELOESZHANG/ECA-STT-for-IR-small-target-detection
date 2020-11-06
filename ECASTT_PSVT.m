function [tenB, tenT,tenN] = ECASTT_PSVT(tenD,lambda1,lambda2,lambda3,mu,wepara,W,C)
%% initialize
rho = 1.5;
tol = 10^(-7);
cg_tol=1e-4;cg_iter=10;
[h,w,t] = size(tenD);
d=2;L=9;B=t;   
r = rankN(tenD,0.1);  %%truncated rank

tenB = tenD;
tenT = zeros(size(tenD));
tenN = zeros(size(tenD));
tenY1 = zeros(size(tenD));tenY3 = zeros(size(tenD));tenY2 = zeros(h,w,d,L,B);
normD = norm(tenD(:));
reweightTen = ones(size(tenD)) .* wepara;
[idx,n,I]=compute_Nsum_indices(C);
bc='symmetric';

iter = 0;
converged = false;
while ~converged
       iter = iter + 1;
       %% Update Z
        tenA = tenB-tenY1/mu;
        tenZ = prox_pstnn(tenA,r,mu);
       %% Update L
        JB = Jacobian(tenB,W,C,bc);
        tenL = softThres(JB-tenY2/mu,lambda1*ones(size(JB))/mu);%
       %% Update B
        JL = AdjJacobian(tenL,W,idx,n,I,bc);
        JY2 = AdjJacobian(tenY2,W,idx,n,I,bc);
        temp1 = @(tenB)JTJ(tenB,W,idx,n,I,bc)+2*tenB;
        temp2 = (tenY1+tenY3+JY2+mu*(tenZ+JL+tenD-tenT-tenN))/mu;
        tenB = cg_solver(temp1,temp2,tenB,cg_iter,cg_tol);
       %% Update T
       tenT = prox_non_neg_l1(tenD-tenB-tenN+tenY3/mu, lambda2*reweightTen/mu);
       %% Update W_ReTE
       reweightTen = 1./ (abs(tenT) + 0.01) .* wepara;
       %% Update N
       tenN = (mu*(tenD-tenB-tenT)+tenY3)/(mu+2*lambda3);
       %% Update Y
        tenY1 = tenY1+mu*(tenZ-tenB);
        tenY2 = tenY2+mu*(tenL-Jacobian(tenB,W,C,bc));
        tenY3 = tenY3+mu*(tenD-tenB-tenT-tenN);
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

function r = rankN(X, ratioN)
    [~,~,n3] = size(X);
  if n3 == 1
    [~, S, ~] = svd(X, 'econ');
    [desS, ~] = sort(diag(S), 'descend');
    ratioVec = desS / desS(1);
    idxArr = find(ratioVec < ratioN);
    if idxArr(1) > 1
        r = idxArr(1) - 1;
    else
        r = 1;
    end        
  else
    D = Unfold(X,n3,1);
    [~, S, ~] = svd(D, 'econ');
    [desS, ~] = sort(diag(S), 'descend');
    ratioVec = desS / desS(1);
    idxArr = find(ratioVec < ratioN);
    if idxArr(1) > 1
        r = idxArr(1) - 1;
    else
        r = 1;
    end
  end