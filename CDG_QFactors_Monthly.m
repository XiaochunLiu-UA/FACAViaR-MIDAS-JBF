%%%%% Estimate Chen, L., Dolado, J., Gonzalo, J. (2021) Quantile Factor Models. Econometrica 89(2):
% 875-910


clear all;
fpath = '...\CDG_QFactors\One Factor\';


[fred_md,name0] =xlsread('FRED_MD_Clean.xlsx');
dat= fred_md(:,2:end);
dates=fred_md(:,1);
datesString=name0(2:end,1);

target=2004+(12-1)/12;
[~, st] = min(abs(dates - target));

totalObs          =size(fred_md,1);
outSampleSize = totalObs - st+1;
tol  = 1e-3;
rmax = 2;
r    = rmax;

TAU= [0.01,0.05,0.1,0.5,0.9,0.95,0.99];


for j = 1:length(TAU)
    tau = TAU(j);   
    parfor i = 1:outSampleSize    
        fprintf('Estimating fators at i = %i \n\n',i); 
        %%%%%%%%%
        %tau = 0.01;
        datt=dat(1:(st+i-1),:);
        [Fhat,Lhat] = IQR(datt,rmax,tol,tau);
        cr=corrcoef(Fhat,datt(:,6));
        if cr(1,2)<0 
            Fhat=-Fhat;
        end
        fname = "lpfactor_" + string(tau) + "_" + i + ".csv";
         fname = fpath + fname;
        writematrix(Fhat, fname);       

         fname1 = "lpfactor_Loadings_" + string(tau) + "_" + i + ".csv";
         fname1 = fpath + fname1;
        writematrix(Lhat, fname1);   

    end
end

  