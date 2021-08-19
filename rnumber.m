%% generate an array of random numbers


days=200;
maxstep=days*2400;
maxtimes=2000;
rng('shuffle');         % set the seed

rnumber11=normrnd(0,1,[1,(maxstep+1)*maxtimes]);

save('rnumber1.mat','rnumber11','-v7.3');

