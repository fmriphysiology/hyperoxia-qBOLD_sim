function figure_dbv_noise

rng('default'); % reset random number generator

OEF=0.4; %assume an OEF of 40%
DBV=0.02; %assume a DBV of 2%

numerator=OEF.*DBV; %i.e. OEF = const x R2' / DBV 

randnum=randn(1000,1).*DBV;

figure;

subplot(1,3,1);
denominator=randnum+DBV; %i.e. DBV + noise (std=DBV)
histogram(numerator./denominator.*100,[0:1:100],'normalization','probability');
xlim([0 100])
ylim([0 0.04])
set(gca,'ytick',[0:0.01:0.04])
set(gca,'xtick',[0:20:100])
axis square;
box on;
grid on;

subplot(1,3,2);
denominator=randnum./5+DBV; %i.e. DBV + noise (std=DBV/5)
histogram(numerator./denominator.*100,[0:1:100],'normalization','probability');
xlim([0 100])
ylim([0 0.07])
set(gca,'ytick',[0:0.02:0.06])
set(gca,'xtick',[0:20:100])
axis square;
box on;
grid on;

subplot(1,3,3);
denominator=randnum./10+DBV; %i.e. DBV + noise (std=DBV/10)
histogram(numerator./denominator.*100,[0:1:100],'normalization','probability');
xlim([0 100])
ylim([0 0.12])
set(gca,'ytick',[0:0.04:0.12])
set(gca,'xtick',[0:20:100])
axis square;
box on;
grid on;
