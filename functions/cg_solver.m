function [x,delta,i]=cg_solver(A,b,x,iter,tol)

% Implements the conjugate gradient method to solve the linear system Ax=b, 
% where A is a function handle for the operation Ax.
%
% Inputs parameters are:
%    b    : particular solution b
%    x    : starting value for x
%    iter : maximum CG iterations
%    tol  : threshold for stoping CG
% Output parameters are:
%    xnew  : solution of conjudate gradient
%    delta : squared norm of the residual r=A(x)-b.
%    i     : number of iterations used


i=0;
r=b-A(x);
d=r;
delta=r(:)'*r(:);

while(i < iter && delta > tol)
  % q=A*d; 
  
  q=A(d);
  alfa=delta/(d(:)'*q(:));
  
  x=x+alfa.*d;
  
  if ~mod(i,50)
    r=b-A(x);
  else
    r=r-alfa*q;
  end
  
  delta_old=delta;
  delta=r(:)'*r(:);
  beta=delta./delta_old;
  d=r+beta.*d;
  i=i+1;
end