function [x,v,s,h,f] = core_bifur

curdir = pwd;
init;
cd(curdir);
opt = contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',500000);
opt=contset(opt,'MinStepsize',1);
opt=contset(opt,'MaxStepsize',100);
opt=contset(opt,'Eigenvalues',1);

% Degradation rate:
ks = 0.125;  ku200 = 0.05;   kmz = 0.5;   kz = 0.1;
kmsl=0.5; ksl=0.1155;
% Transcription rate:
gs = 18000; gu200 = 2100;   gmz = 11;   gz = 100; 
gmsl=90; gsl=50000;
% Hills function threshold :
 I0s=100000; z0u200 = 220000;   z0mz = 25000;   s0u200 = 180000;   s0mz = 180000; u2000 = 10000;
sl0u200=220000; sl0s=225000; s0msl=180000; s0s=300000;
% Cooperativity:
nsmz = 2;  nIs = 2; nzu200 = 3;   nsu200 = 2;   nzmz = 2;   nu200 = 6; 
nslu200=1;nsls=3; nsmsl=1; nss=5;
% fold change
lamdazu200 =0.1;   lamdasu200 = 0.1;  lamdazmz = 7.5;   lamdasmz = 10; lamdaIs=3;
lamdaslu200=0.4; lamdasls=0.5; lamdasmsl=0.5; lamdass=0.4;
% external signal
I=0;

ap = 1; %describes the index of parameter for which the bifurcation is drawn using the init_EP_EP function. Currently, ap=1, thus bifurcation parameter is s (SNAIL levels)
handles = feval(@core);
tspan = 0:100:50000;

% initial condition
x_start = [33554.833280 56.5 0 0 0 0];

%calculating steady state for given initial condition 
[t,x_time] = ode15s(@(t,kmrgd)handles{2}(t,kmrgd,I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl),tspan,x_start);
x_init = x_time(end,:)';

%drawing bifurcation using a continuation method
[x0,v0] = init_EP_EP(@core,x_init,[I,lamdazu200,nzu200,z0u200,nu200,u2000,lamdasu200,nsu200,s0u200,lamdazmz,nzmz,z0mz,lamdasmz,nsmz,s0mz,lamdaslu200,sl0u200,nslu200,lamdasmsl,nsmsl,s0msl,lamdasls,nsls,sl0s,lamdass,nss,s0s,lamdaIs,I0s,nIs,gu200,ku200,gmz,kmz,gz,kz,gmsl,kmsl,gs,ks,gsl,ksl],ap);
[x,v,s,h,f] = cont(@equilibrium, x0, v0,opt);