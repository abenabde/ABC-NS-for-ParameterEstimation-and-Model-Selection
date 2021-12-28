clc;clear all;close all;
%% Matlab Script for Post Processing
% Plot the results 
% M_1 : Linear model composed of two parameters
% M_2 : Cubic model composed of three parameters
load Output_ABCNS_34
% Histograms of the parameters of the selected model 
% The selected model is the cubic model
pop = population_2(:,3*sim-2:3*sim) ;  % extract particles from the last population
true = [0.05 20 1000]; % true values
gra = colormap();
f = {'\fontname{Times} \fontsize{12} \itc' '\fontname{Times} \fontsize{12} \it{k}' ...
    '\fontname{Times} \fontsize{12} \it{k_3}'};
last = pop; % the last population is stored in pop
yu = [1 2 3.5];
mean_value = mean(last);
gra = colormap(gray);
figure(1)
for km = 1 : 3
    subplot(2,2,yu(km));
    hist(pop(:,km));
    hold on
    plot(mean_value(1,km),0,'g^','MarkerSize',5,'MarkerFaceColor','g',...
        'MarkerEdgeColor','g')
    set(gcf,'Color','w')
    h = findobj(gca,'Type','patch');
    h.FaceColor = [gra(50,:)];
    h.EdgeColor = 'k';
    xlabel((f(1,km)))
    ylabel('\fontname{Times} \fontsize{12} Frequency')
end
%----------------------------------------------------------
% This figure shows the posterior model probabilities over some populations
proba = [ix(1:sim,:)./1000];
proba = [0.5 0.5;proba];
is = [(1:3:length(proba))']; % select some populations
figure(2)
for j=1:length(is)
subplot(4,3,j)
bar(proba(is(j),1:end),'FaceColor',[gra(50,:)],'EdgeColor','k')
axis([0 3 0 max(proba(is(j),1:end))+0.1])
set(gcf,'Color','w')
title(['$p$='  num2str(is(j)) ,',' ,'$~\varepsilon$=' num2str(threshold_vec(is(j)),'%2.2e\n')],'FontSize',9,'FontName','Times','Interpreter','latex')
x_tick_name={'$$\mathcal{M}_1$$','$$\mathcal{M}_2$$'};
set(gca,'TickLabelInterpreter','lATEX')
set(gca,'XTickLabel',x_tick_name);
end
%-------------------------------------------------
um_1  =  c_simulate(1,0.05,50,1e3); % Cubic model (c,k,k3)
um_pred  =  c_simulate(1,mean_value(1,1),mean_value(1,2),mean_value(1,3)); 
figure(3)
plot(um_1,'k-','LineWidth',2)
hold on
plot(um_pred,'b-.','LineWidth',1.5)
xlabel('\fontname{Times} \fontsize{12} Discrete time points')
ylabel('\fontname{Times} \fontsize{12} Displacement')
axis([0 500 -0.2 0.2])
legend('\fontname{Times} \fontsize{13} Observed data',...
    '\fontname{Times} \fontsize{13} Predicted data',...
     'Location','Best')