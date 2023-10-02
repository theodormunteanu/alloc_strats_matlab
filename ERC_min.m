function [str,vol,RC] = ERC_min(sigs,corr)
   % ERC portfolio when colatilities and correlation matrix is given
   % 
   cov_mat = (sigs'.*sigs).*corr;
   n = length(sigs);
   func = @(x) sqrt(x*cov_mat*x')-1/n*sum(log(x));
   str = fmincon(func,ones(1,n)/n,[],[],ones(1,n),1,zeros(1,n),ones(1,n));
   vol = sqrt(str*cov_mat*str');
   RC = str.*(cov_mat*str')'/vol;
end