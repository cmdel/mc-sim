clf;
paths=50;
dt=0.003937; %254 days per year
lambda = 5;
years=1;
xi=0.6;
S0=5;
mu=1;
days=floor(years/dt);

subplot(2,1,1)
set(gcf, 'Position', get(0,'Screensize'))
CIRpaths(S0,mu,lambda,xi,years,dt,paths);
xlabel('Time(business days)');
ylabel('S(t)');
xlim([0 days]);
ylim([-2 6]);
title(gca,['Euler Discritisation, sampling from N_{(0,1)} for ',num2str(years),' year(s)']);
hk=line([0 days],[1 1],'Color','k');
set(hk,'LineStyle','--');
hr=line([0 days],[0 0],'Color','r');
set(hr,'LineStyle','--');
    
subplot(2,1,2)
hold all
for t=1:paths
    plot(cirpath(0:dt:years,lambda,mu,xi,S0),'LineWidth',2);
end
xlabel('Time(business days)');
ylabel('S(t)');
xlim([0 days]);
ylim([-2 6]);
title(gca,['Euler Discritisation, sampling from the terminal non-central \chi^2 distribution for ',num2str(years),' year(s)']);
hk=line([0 days],[1 1],'Color','k');
set(hk,'LineStyle','--');
hr=line([0 days],[0 0],'Color','r');
set(hr,'LineStyle','--');

fprintf('Completed\n');