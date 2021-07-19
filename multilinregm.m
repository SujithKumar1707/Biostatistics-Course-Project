tifunction out = multilinregm(indat, alpha)
out= struct('alpha',0,'betacap',0,'reg_dof',0,'err_dof',0,'tot_dof',0,'nparams',0,'SE',0,'param_range',0,'xmat',0,'ymat',0,'SSE',0,'SSR',0,'SST',0,'Fr',0,'RMSE',0,'Rsq',0,'Rsqadj',0,'sigmacapsq',0,'H0_F',0,'tcals',0,'H0_beta',0,'tpvalues',0,'fpvalue',0);
n=length(indat(:,1)); nx=length(indat(1,:))-1; out.nparams=nx+1; out.reg_dof=nx;
out.xmat=zeros(n,nx+1); out.ymat=zeros(n,1);out.err_dof=n -out.nparams; out.tot_dof=n-1;
out.alpha=alpha;
for i=1:1:n
    out.xmat(i,1) =1;
    for k=2:1:nx+1
        out.xmat(i,k)=indat(i, k-1);
    end
    out.ymat(i,1)=indat(i,nx+1);
   
end
out.betacap=((out.xmat'*out.xmat)^-1)*(out.xmat'*out.ymat);
ycap=out.xmat*out.betacap;
out.SSE=sum((out.ymat-ycap).^2);
out.sigmacapsq=out.SSE/out.err_dof;
out.covbeta=((out.xmat'*out.xmat)^-1)*out.sigmacapsq;
out.SST=out.ymat'*out.ymat-sum(out.ymat)^2/n; out.SSR=out.SST-out.SSE;
out.Fr = (out.SSR/out.reg_dof)/(out.SSE/out.err_dof);
out.fpvalue=fcdf(out.Fr,out.reg_dof,out.err_dof,'upper');
if out.fpvalue>1, out.fpvalue=out.fpvalue/2; end
Ftablow=finv(alpha/2,out.reg_dof,out.err_dof); Ftabhigh=finv(1-alpha/2, out.reg_dof,out.err_dof);
if out.Fr>Ftablow && out.Fr<Ftabhigh
    out.H0_F=0;
else
    out.H0_F=1;
end
out.tcals=zeros(nx+1,1); out.tpvalues=zeros(nx+1,1); out.SE=zeros(nx+1,1);
for i=1:1:nx+1
    out.tcals(i,1)=out.betacap(i)/sqrt(out.covbeta(i,i));
     out.SE(i,1)=sqrt(out.covbeta(i, i));
     out.tpvalues(i,1)=2*tcdf(abs(out.tcals(i,1)), out.err_dof, 'upper');
     if out.tpvalues(i,1)> 1, out.tpvalues(i,1)=out.tpvalues(i,1)/2;end
end
out.H0_beta=zeros(nx+1,1);
ttablow=tinv(alpha/2,out.err_dof); ttabhigh=tinv(1-alpha/2, out.err_dof);
for i=1:1:nx+1
    if out.tcals(i)>ttablow && out.tcals(i)<ttabhigh 
        out.H0_beta(i) = 0;
    else 
        out.H0_beta(i) = 1;
    end
end
out.param_range=zeros(nx+1,1);
for i=1:1:nx+1
    out.param_range(i,1)=abs(tinv(alpha/2,out.err_dof))*sqrt(out.covbeta(i,i));
    %out.param_range(i,2)=out.betacap(i)+abs(tinv(alpha/2,out.err_dof))*sqrt(out.covbeta(i,i));
end
out.Rsq = out.SSR/out.SST;
out.Rsqadj=1-(out.SSE/out.err_dof)/(out.SST/out.tot_dof); 
out.RMSE=sqrt(out.sigmacapsq);
end


     


    


