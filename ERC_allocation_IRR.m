function [str,sig] = ERC_allocation_IRR(cov_mat,durations,convexities,varargin)
   % we compute the ERC allocation when only interest rate risk is being
   % considered, and bonds durations and convexities are computed. 
   
   p = inputParser;
   addOptional(p,'CI',[],@(x) isvector(x));
   addOptional(p,'Reduction',[],@(x) isnumeric(x));
   addOptional(p,'bench',[],@(x) isvector(x));
   
   n = length(durations);
  %%
   if isempty(convexities)==1
      cov_mat_tilde = (durations'*durations).*cov_mat;
   else
      cov_mat_tilde = (durations'*durations).*cov_mat + 1/2*(convexities'*convexities).*(cov_mat).^2;
   end
   vol = @(x) sqrt(x*cov_mat_tilde*x');
   RC = @(x) x.*(cov_mat_tilde*x')'/vol(x);
   str = fmincon(@(x) norm(RC(x)-ones(1,n)/n,2),ones(1,n),[],[],ones(1,n),1);
   sig = vol(str);
end

