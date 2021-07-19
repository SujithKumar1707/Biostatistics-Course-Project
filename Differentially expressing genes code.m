
% mx=1;
% my=200;
% vx=200;
% vy=2;
alpha= 0.01;
% dax= mx + sqrt(vx) * randn(100,1);
% day= my+ sqrt(vy) * randn(100,1);
dax1=[24.1377,6.99386,3.51592,-6.57815,5.74084,12.4085,0.286943,19.4414,5.14252,22.4605];
day1=[67.9946,60.9062,58.2787,72.5177,61.1863,52.2541,71.3123,77.256,68.8314,68.0983];
dax2=[48.3578,52.6324,39.7174,62.0014,60.1729,51.4901,44.3322,51.0155,58.6339,78.0409];
day2=[22.2971,18.5077,36.6209,42.3734,31.9061,28.7108,26.0777,14.5733,25.1358,26.8633];
dax3=[17.1639,13.3421,23.3693,12.2059,18.1567,3.19916,32.3981,0.765049,3.84458,6.53739];
day3=[50.8454,48.1654,53.3688,40.5515,44.1172,35.5211,51.6407,48.8935,44.2119,36.5952];
dax4=[11.8076,9.64813,0.444355,18.0779,12.475,22.2247,18.6849,25.9072,12.1507,14.0054];
day4=[124.215,103.338,82.5999,90.5018,102.036,109.109,98.3123,116.977,72.302,110.605];
dax5=[64.3253,52.3742,35.5534,46.1861,65.8666,24.0491,27.3642,39.9469,36.8725,35.7628];
day5=[125.972,104.277,104.114,60.4773,110.236,70.2754,100.286,104.291,108.94,106.84]; 
%out=ttest(dax3,day3,alpha);
outm=general_ttestmat(dax5,day5,alpha);
% [a,b,c,d] = vartest2(dax,day);
% [a,b,c,d] = ttest2(dax,day)
%disp(out);
disp(outm);
function out= general_ttestmat(x,y,alpha)
out=struct('mean_x', 0,'mean_y',0,'var_x',0,'var_y',0, 'dof_x', 0, 'dof_y',0,'tot_dof',0,'F_cal',0 ,'var_r_l',0,'var_r_h',0,'diff_mn_l',0,'diff_mn_h',0,'sd',0,'tcal', 0, 'nh_F',0 ,'nh_t',0,'fpval',0,'pvalt', 0,'tpval',0);
out.mean_x = mean(x);
out.mean_y = mean(y);
out.var_x = var(x);
out.var_y = var(y);
out.dof_x = length(x)-1;
out.dof_y = length(y)-1;
[a,b,c,d] = vartest2(x,y,'Alpha',alpha);
out.nh_F=a;
out.fpval = b;
out.var_r_l= c(1); 
out.var_r_h = c(2) ; 
out.F_cal=d.fstat;
if out.nh_F == 0
   [a,b,c,d] = ttest2(x,y,'Alpha',alpha,'Vartype','equal');
   out.nh_t=a;
   out.tpval = b;
   out.diff_mn_l= c(1); out.diff_mn_h = c(2) ; out.tcal=d.tstat;
   out.tot_dof = d.df;
   out.sd = d.sd;
else 
   [a,b,c,d] = ttest2(x,y,'Alpha',alpha,'Vartype','unequal');
   out.nh_t=a;
   out.tpval = b;
   out.diff_mn_l= c(1); out.diff_mn_h = c(2) ; out.tcal=d.tstat;
   out.tot_dof = d.df;
   out.sd = d.sd;
end
end