%% script for plotting data in simulationData.mat
clc; clear;
load('simulationData.mat');
set(0,'defaultLineLineWidth',1,'DefaultLineMarkerSize',5)

%% Plot the trajectories for selected simulations

plot_indices = {[1 7 16]
                [3 5 7 8]
                [9:3:26]};
titles = {'','DAMPC simulations','PE simulations'};
names = {'AllTraj','DAMPC','PE'};
colors = {'b','r','k','g','m','y'};
markers = {'o','s','*','^','v','>'};
legends = {{'PAMPC','DAMPC-best','PE-best'},
            {'$\hat{N}=2$','$\hat{N}=4$','$\hat{N}=6$','$\hat{N}=7$'},
            {'$\rho=0.01$','$\rho=0.07$','$\rho=0.2$','$\rho=0.5$','$\rho=1$','$\rho=4$'}}

Tsim = length(data(1).x)-1;
x_ref = data(1).x_ref;

f = 80;
for plots = 1:3
h = figure(f); f = f+1;
clf;
i = 1;
for ind=plot_indices{plots}
    subplot(4,2,1)
    title(titles{plots})
    hold on;
    plot(0:Tsim,x_ref(1,:),'k--')
    plot(0:Tsim,data(ind).x(1,:),['-',markers{i},colors{i}])
    xlim([0 Tsim])
    ylim([-1 3])
    ylabel('$x_1$')

    subplot(4,2,3)
    hold on;
    plot(0:Tsim,x_ref(2,:),'k--')
    plot(0:Tsim,data(ind).x(2,:),['-',markers{i},colors{i}])
    xlim([0 Tsim])
    ylim([-3 1])
    ylabel('$x_2$')

    subplot(4,2,5)
    hold on;
    plot(0:Tsim-1,data(ind).u(1,:),['-',markers{i},colors{i}])
    xlim([0 Tsim])
    ylim([-4 4])
    ylabel('$u_1$')
    legend(legends{plots});
    
    subplot(4,2,7)
    hold on;
    plot(0:Tsim-1,data(ind).u(2,:),['-',markers{i},colors{i}])
    xlim([0 Tsim])
    ylim([-4 4])
    ylabel('$u_2$')

    subplot(4,2,[2,4,6,8]); 
    hold on; 
    plotregionLine(-data(ind).cont.H_theta,-data(ind).cont.h_theta_k,[],[],colors{i},[],[],[],[]);
    plot(0.95,0.3,'k*')
    legend({'Bounds','$\theta^*$'})
    xlabel('$\theta_1$')
    ylabel('$\theta_2$')
    xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
 
 i = i+1;
end
saveas(h,names{plots},'epsc')
end

%% Plot J-rho for the PE part

i = 1;
% extract rho_PE
for ind = indices.PE
    rho_PE(i) = data(ind).cont.rho_PE;
    i = i+1;
end
% extract J
J = [data.J];

h = figure(f); f = f+1;
clf;
subplot(1,2,1)
semilogx(rho_PE,J(indices.PE),'r*','MarkerSize',8);
xlabel('$\rho$')
ylabel('$J$')
ylim([50,400]);
title('PE + PAMPC')
subplot(1,2,2)
plot([0:7],J([1 indices.DAMPC]),'b*','MarkerSize',8);
xlabel('$\hat{N}$')
ylabel('$J$')
ylim([50,400]);
title('DAMPC')
saveas(h,'Comparison','epsc')
%% Animate evolution of the parameter set
% Choose indices
vid_indices = [1 5 7 12 15 24];
vid_names = {'PAMPC','DAMPC_4','DAMPC_6','PE_7e2','PE_2e1','PE_4'}
i = 2;
j = 1;
for ind = vid_indices
    x = NaN*ones(2,Tsim+1);
    u = NaN*ones(2,Tsim);
    x(:,1) = data(ind).x(:,1);       
    h = figure(f); f = f+1;
    for k = 1:Tsim
        x(:,k+1) = data(ind).x(:,k+1);
        u(:,k) = data(ind).x(:,k);
        
        clf;
        subplot(4,2,1)
    %     title(titles{plots})
        hold on;
        plot(0:Tsim,x_ref(1,:),'k--')
        plot(0:Tsim,x(1,:),['-',markers{i},colors{i}])
        xlim([0 Tsim])
        ylim([-1 3])
        ylabel('$x_1$')

        subplot(4,2,3)
        hold on;
        plot(0:Tsim,x_ref(2,:),'k--')
        plot(0:Tsim,x(2,:),['-',markers{i},colors{i}])
        xlim([0 Tsim])
        ylim([-3 1])
        ylabel('$x_2$')

        subplot(4,2,5)
        hold on;
        plot(0:Tsim-1,u(1,:),['-',markers{i},colors{i}])
        xlim([0 Tsim])
        ylim([-4 4])
        ylabel('$u_1$')
    %     legend(legends{plots});

        subplot(4,2,7)
        hold on;
        plot(0:Tsim-1,u(2,:),['-',markers{i},colors{i}])
        xlim([0 Tsim])
        ylim([-4 4])
        ylabel('$u_2$')

        subplot(4,2,[2,4,6,8]); 
        hold on; 
        plotregionLine(-data(ind).cont.H_theta,-data(ind).h_theta_k_mat(:,k),[],[],colors{i},[],[],[],[]);
        plot(0.95,0.3,'k*')
        legend({'Bounds','$\theta^*$'})
        xlabel('$\theta_1$')
        ylabel('$\theta_2$')
        xlim([-1 1]);ylim([-1 1]);zlim([-1 1]);
        
        movieFrames(k) = getframe(h);
    end
    vidName = ['Simulation',num2str(ind)];
    myWriter = VideoWriter(vid_names{j});
    myWriter.FrameRate = 2;
    open(myWriter);
    writeVideo(myWriter, movieFrames);
    close(myWriter);
j = j+1;
end

%% Runtimes

% extract rho_PE
for ind = 1:length(data)
    runtimes(:,ind) = data(ind).runtimes';
end
mean_run = mean(runtimes,1);