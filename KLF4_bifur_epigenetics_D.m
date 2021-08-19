clear all
[x,v,s,h,f] = KLF4_bifur_epigenetics; 
a = x(11,:); %bifurcation parameter
b = x(2,:); 
c = a./1000;

%% Based on eigenvalues to judge stable vs. unstable states
ind = zeros(1,4);
snum = size(f);
num = snum(2);
j = 1;
n = 1;

for n = 1:1:(num-1)
    x1 = find(f(:,n) > 0);
    x2 = find(f(:,n+1) > 0);
    if isempty(x1) && ~isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    elseif ~isempty(x1) && isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    end
end

%%

amat=[c(ind(1)) c(ind(2)) c(ind(5)) c(ind(end))];
figure1 = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
% % 
% % ax2 = subplot(2,2,2);
plot(c(1:ind(1)),b(1:ind(1)),'b');
hold on
plot(c(ind(1)+1:ind(2)),b(ind(1)+1:ind(2)),'r');
plot(c(ind(2)+1:ind(7)),b(ind(2)+1:ind(7)),'b'); 
plot(c(ind(7)+1:ind(end)),b(ind(7)+1:ind(end)),'r');
plot(c(ind(end)+1:end),b(ind(end)+1:end),'b');
% 
% 
xlim([0 200]);
% % xlabel('I ext (10^3 molecules)');
% ylabel('ZEB mRNA');
% sound(sin(1:3000));
