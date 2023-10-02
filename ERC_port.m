function [str,vol,RC] = ERC_port(sigs,corr,varargin)

 % Find the ERC (Equally Risk Charged) portfolio given carbon reduction
 % portfolio. 
  p = inputParser;
  addOptional(p,'CI',[],@(x) isvector(x));
  addOptional(p,'Reduction',[],@(x) isnumeric(x));
  addOptional(p,'bench',[],@(x) isvector(x));
  parse(p,sigs,corr,varargin{:});
  cov_mat = (sigs'.*sigs).*corr;
  n = length(sigs);
  f = @(x) sum(x.^2 .* (cov_mat*x')'.^2);
  g = @(x) x*cov_mat*x';
  
  CIs = p.Results.CI;R=  p.Results.R;bench = p.Results.bench;
  if isempty(CIs)==1 && isempty(R)==1 && isempty(bench)==1
     str = fmincon(@(x) f(x)/g(x)^2,ones(1,n)/n,[],[],ones(1,n),1,zeros(1,n),ones(1,n));
  else
     C = CIs;D = (1-R)*CIs*bench';
     str = fmincon(@(x) f(x)/g(x)^2,ones(1,n)/n,C,D,ones(1,n),1,zeros(1,n),ones(1,n));
  end
  vol = sqrt(str*cov_mat*str');
  RC = str.*(cov_mat*str')'/vol;
end