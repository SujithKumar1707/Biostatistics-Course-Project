clc
format long;
ddataset = xlsread('CH19B012.csv');
X = ddataset(:,1);
Y = ddataset(:,2);

initx = X(1);
inity = Y(1);

delx = X(2) - X(1);
maxx = X(length(X(:,1)),1);

lambda = 50;
lambda_s = lambda;
lambda_l = 1e-3;
lambda_h = 1e8;
tune_p = 2;
nite = 150000;
stat = 1;
pas = zeros(nite,5);

dof = length(X)-5;
alpha = 0.02563;
delp = 1e-6;
tol = 1e-12;
param = zeros(5,1);
param(1,1) = 1.0;
param(2,1) = -1.7;
param(3,1) = 0.5;
param(4,1) = 1.0;
param(5,1) = 0.2;

out = struct('alpha',alpha,'params',0,'total_iterations',0,'initial_lambda',0,'RMSE',0,'final_lambda',0,'tune_param',0,'tpvalues',0,'errDOF',0,'totDOF',0,'regDOF',0,'MSE',0,'nparam',0,'tstat',0,'covmat',0,'SE_params',0,'H0_params',0,'SST',0,'SSR',0,'SSE',0,'fcal',0,'H0_f',0,'fpvalue',0,'R_square',0,'R_square_adj',0);

out.errDOF = dof;
out.initial_lambda = lambda;
out.totDOF = length(X)-1;
out.regDOF = 4;
out.tune_param = tune_p;
out.nparam = 5;

diffa = @(a,b,c,d,h,delx,maxx,delp) ((expy(a+delp,b,c,d,h,delx,maxx)-expy(a,b,c,d,h,delx,maxx))/delp);
diffb = @(a,b,c,d,h,delx,maxx,delp) ((expy(a,b+delp,c,d,h,delx,maxx)-expy(a,b,c,d,h,delx,maxx))/delp);
diffc = @(a,b,c,d,h,delx,maxx,delp) ((expy(a,b,c+delp,d,h,delx,maxx)-expy(a,b,c,d,h,delx,maxx))/delp);
diffd = @(a,b,c,d,h,delx,maxx,delp) ((expy(a,b,c,d+delp,h,delx,maxx)-expy(a,b,c,d,h,delx,maxx))/delp);
diffh = @(a,b,c,d,h,delx,maxx,delp) ((expy(a,b,c,d,h+delp,delx,maxx)-expy(a,b,c,d,h,delx,maxx))/delp);

for w=1:nite
    if stat ==1
        jmat = [diffa(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx,delp) diffb(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx,delp) diffc(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx,delp) diffd(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx,delp) diffh(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx,delp)];
        jtjm = jmat'*jmat;
        dy = Y - expy(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx);
        if lambda_s ~= 0 
            if w ==1
                chia = sum(dy.^2);
            end
        else
            chia = sum(dy.^2);
        end
    end
    pas(w,:) = param(:,1);
    if lambda_s == 0
        for i = 1:length(jtjm),jtjm(i,i) = jtjm(i,i) + 1e-1;end
    else
        for i = 1:length(jtjm),jtjm(i,i) = jtjm(i,i) + jtjm(i,i)*lambda;end
    end
    ijtjm = jtjm^-1;
    
    dyn = Y - expy(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx);
    nparam = param + ijtjm * (jmat' *dyn);
    
    dym = Y - expy(nparam(1,1),nparam(2,1),nparam(3,1),nparam(4,1),nparam(5,1),delx,maxx);
    chib = sum(dym.^2);
    
    if abs(chia-chib) <tol || w==nite
        param = nparam;
        SSE = chib;
        SST = sum(Y.^2) - sum(Y)^2/length(Y);
        SSEr = SSE/dof;
        RMSE = sqrt(SSEr);
        out.RMSE = RMSE;
        SSTr = SST/(length(Y)-1);
        rsquare = 1 - SSE/SST;
        rsqadj = 1-((length(Y)-1)/dof)*(1-rsquare);
        cov_mat = ijtjm*SSEr;
        param_err = zeros(5,2);
        param_err(:,1) = param - (abs(tinv(alpha/2,dof)) * sqrt(diag(cov_mat)));
        param_err(:,2) = param + (abs(tinv(alpha/2,dof)) * sqrt(diag(cov_mat)));
        tmp = sqrt(diag(cov_mat));
        
        out.tpvalues = zeros(5,1);
        for ip = 1:5
            out.tpvalues(ip,1) = tpval(param(ip)/tmp(ip), out.errDOF);
        end
        out.SE_params = param_err;
%         disp(['param (a,b,c,d,h), ' 'param_low, ' 'param_high ']);
        disp([param param_err]);
        disp(['total iterations, ' 'SSE, ' 'RMSE, ' 'rsquare, ' 'rsqadj, ' 'lambda']);
        out.total_iterations = w;
        disp([w SSE RMSE rsquare rsqadj lambda]);
        out.params = param;
        out.SSE = SSE; out.SST = SST; out.covmat = cov_mat;
        out.SSR = out.SST-out.SSE;
        out.final_lambda = lambda;
        out.R_square = rsquare;
        out.MSE = SSEr;
        out.R_square_adj = rsqadj;
        out.fcal = (out.SSR/out.regDOF)/SSEr;
        out.fpvalue = fpval(out.fcal,out.regDOF,out.errDOF);
        
        if out.fpvalue < alpha
            out.H0_f = 1;
        else
            out.H0_f = 0;
        end
        
        out.tstat = param.*sqrt(diag(cov_mat)).^-1;
        out.tpvalues = zeros(length(out.tstat),1);
        out.H0_params = zeros(length(out.tstat),1);
        for q=1:length(out.tstat)
            out.tpvalues(q) = tpval(out.tstat(q),out.errDOF);
            if out.tpvalues < alpha
                 out.H0_params(q) = 1;
            end
        end
        break;
    end
    if lambda_s ~= 0
        if chib < chia
            if lambda/tune_p >= lambda_l
                lambda = lambda/tune_p; 
            end
            param = nparam;
            chia = chib;
            stat = 1;
        else
            if lambda * tune_p <= lambda_h
                lambda = lambda * tune_p; 
            end
            stat = 0;
        end
    else
        param = nparam;
    end
end

% disp(out);
param_low = param_err(:,1);
std_error = param(:)-param_low(:);


final_answers = table(out.params,std_error,out.tpvalues,'VariableNames',{'Estimate' 'Std Error' 'pvalues'});
disp(final_answers);

disp('Total iteration');
disp(out.total_iterations);
disp('Final SSE');
disp(out.SSE);
disp('Final lambda');
disp(out.final_lambda);
disp('Fcal');
disp(out.fcal);
disp('Covariance value at (3,2)');
disp(out.covmat(3,2));

plot(X,Y,'or'); hold on;
vp = expy(param(1,1),param(2,1),param(3,1),param(4,1),param(5,1),delx,maxx);
plot(X,vp,'-b'); hold off;
cx = gca; cx.FontSize = 16;
xlabel('X'); ylabel('Y');
legend('dataset','nonlinear fit');
disp(length(X));


function pval = fpval(f, adof, bdof)
    if f < 1
       pval = fcdf(f,adof,bdof);
    else
        pval = fcdf(1/f,bdof,adof);
    end
end

function rpval = tpval(t, v)
rpval = betainc(v/(v+t^2),v/2,0.5);
end