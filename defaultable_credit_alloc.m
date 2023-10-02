function [str,ret,EL] = defaultable_credit_alloc(cov_mat_delta_spread,spreads,SDs,PDs,LGDs,...
           expiries,mkt_prices,FVs,coupons,frequencies,varargin)
    % INPUTS:
    % ------------------------
    % cov_mat_delta_spread = covariance matrix between the changes in
    % credit spread (Required)
    % spreads = vector of credit spreads;
    % SDs = spread durations
    % PDs = annual(ized) probability of defaults 
    % LGDs = Loss given defaults
    
    % OUTPUTS:
    % str = structure of the credit portfolio
    % ret = expected return of the credit portfolio
    % EL = expected loss of the credit portfolio
    cov_mat_tilde = (SDs'*SDs).*cov_mat_delta_spread;
    n = size(cov_mat_tilde,1);
    p = inputParser;
    addOptional(p,'t',1/2);
    addOptional(p,'exp_chgs',zeros(1,n));
    addOptional(p,'IRR_chgs_cov_mat',zeros(n,n));
    parse(p,varargin{:});
    t = p.Results.t;
    
    exp_chgs = p.Results.exp_chgs;
    IRR_chgs_cov_mat = p.Results.IRR_chgs_cov_mat;
    if IRR_chgs_cov_mat == zeros(n,n)
        % interest rate risk is hedged
        func = @(x)sqrt(x*cov_mat_tilde*x');
        str = fmincon(func,ones(1,n)/n,[],[],ones(1,n),1,zeros(1,n),ones(1,n));
        ret = str*spreads'*t-(str.*spreads)*exp_chgs'-t*(PDs.*LGDs)*str';
    else
        % Interest rate risk is not hedged
        bond_yields = arrayfun(@(i) bond_yield(FVs(i),coupons(i),frequencies(i),mkt_prices(i),expiries(i)),...
            1:length(SDs));
        for i = 1:length(SDs)
            [~,EDs(i)] = bond_price(FVs(i),coupons(i),frequencies(i),bond_yields(i),expiries(i));
        end
        cov_mat_tilde = (SDs'*SDs).*cov_mat_delta_spread+(EDs'*EDs).*IRR_chgs_cov_mat; 
        func = @(x) sqrt(x*cov_mat_tilde*x');
        str = fmincon(func,ones(1,n)/n,[],[],ones(1,n),1,zeros(1,n),ones(1,n));
        ret = str*spreads'*t-str.*spreads*exp_chgs'-t*(PDs.*LGDs)*str';
        
        
    end
    PDs_applied = 1-(1-PDs).^t;
    EL = PDs_applied.*LGDs*str';% expressed in % of invested capital
end

