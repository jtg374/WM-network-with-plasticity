clc;clear all;close all;    
%% load parameters
param = NDF_with_Plasticity_Parameters()


%% unpack Connectivity profile 
MEE = param.MEE;
MEI = param.MEI;
MIE = param.MIE;
MII = param.MII;

%% output time resolution
dt_store = param.dt_store;

%% pack initial values
nx = param.N;
y0 = [0;              % Stimlus Current Strength
      0;              % Wipe Current Strength
      zeros(6*nx,1);  % 6*N state variables
      reshape(MEE,nx*nx,1) % E to E Connection Strength
      ]; 
%% load timing
TStimOn   = param.TStimOn;
TStimOff  = param.TStimOff;
TDelayOff = param.TDelayOff;
Tmax = param.Tmax;

%% Solving ODE equations
clear textprogressbar % will cause trouble if integrate with @odetpbar while textprogressbar is in environment
options = odeset('RelTol',1e-3,'AbsTol',1e-5,'OutputFcn',@odetpbar); % will print progressbar
disp(['Integration started at: ',datestr(now,'HH:MM:SS')])
[t,y] = ode23(@(t,y0) NDF_with_Plasticity_Equations(t,y0,param),...
    0:dt_store:Tmax,y0,options);
disp(['Integration ended at:   ',datestr(now,'HH:MM:SS')])

nt = length(t);
Mt = y(:,nx*6+3:end);MEEt = reshape(Mt,nt,nx,nx);
MEEt = MEEt(TDelayOff/dt_store,:,:);
MEEt = permute(MEEt,[2 3 1]); % put time on 3rd dimention
Rt = y(:,3:nx*6+2);Rt = reshape(Rt,nt,nx,6);
RE = Rt(:,:,1)';RI = Rt(:,:,2)';SEE = Rt(:,:,3)';SIE = Rt(:,:,4)';SEI = Rt(:,:,5)';SII = Rt(:,:,6)'; % put time on 2nd dimention
clear Rt;
% Input_forget=y(:,2);
% Input_stim = y(:,1); 
% Input = Input_stim - Input_forget;
clear y


%% Figures
close all


datapath = ['../../data/FR_Curr_ring_RK4_distractor_with_Plasticity/' datestr(now,'yymmdd_HH_MM_')];
datapath = ['../../data/FR_Curr_ring_RK4_distractor_with_Plasticity/' datestr(now,'yymmdd_HH_MM_')];
mkdir(datapath)

% 

% h2=figure(2); %imagesc([RE RE1])
% tIndex = t>=TStimOn(1) & t<TDelayOff(10);
% subplot(2,1,1);imagesc(RE(:,tIndex));title('first 10 trials')
% ylabel('position (80\theta / 2\pi)','FontSize',10)
% tIndex = t>=param.TStimOn(end-9) & t<TDelayOff(end);
% subplot(2,1,2);imagesc(RE(:,tIndex));title('last 10 trials')
% ylabel('position (80\theta / 2\pi)','FontSize',10)
% xlabel('Time (a.u.)','FontSize',14)
% saveas(h2,[datapath,'/2.fig'])
% saveas(h2,[datapath,'/2.jpg'])

% x = param.x;
% h7=figure(7);hold on;
% firstTestDelay = t>=TStimOff(param.nTrial+1) & t<TDelayOff(param.nTrial+1) ;
% colors = copper(sum(firstTestDelay));
% set(gca,'ColorOrder',colors)
% plot(x,RE(:,firstTestDelay)','LineWidth',0.5)
% title('Evolution of activity pateern at delay of last trial')
% xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% xlim([-pi pi]);
% set(gca,'Ytick',0:50:100,'FontSize',14)
% ylim([0 120])
% hold off
% saveas(h7,[datapath,'/7.fig'])
% saveas(h7,[datapath,'/7.jpg'])


% figure(3);hold on
% i = round(1*pi/dx)+1;
% plot(0:dt_store:Tmax,RE(i,:),'Color',color,'LineWidth',1)
% title('cell at \theta = 0')
% xlabel('Time (s)');ylabel('Firing Rate (Hz)','FontSize',14)
% set(gca,'Xtick',Tinit(21:20:length(Tinit)),'FontSize',14)
% set(gca,'XTickLabel',(20:20:length(Tinit)-1)*(Tstim+Tmemory+Tforget)/1000)
% xlim([Tinit(1)-500 Tmax-500])
% 
% set(gca,'Ytick',0:50:100,'FontSize',14)
% ylim([0 120])
% 
% % firing rate of cell at \theta = \pi
% figure(2);hold on
% i = 1;
% plot(0:dt_store:Tmax,RE(i,:),'Color',color,'LineWidth',1)
% title('cell at \theta = \pi')
% xlabel('Time (s)');ylabel('Firing Rate (Hz)','FontSize',14)
% set(gca,'Xtick',Tinit(21:20:length(Tinit)),'FontSize',14)
% set(gca,'XTickLabel',(20:20:length(Tinit)-1)*(Tstim+Tmemory+Tforget)/1000)
% set(gca,'Ytick',-10:10:40,'FontSize',14)
% ylim([0 25])
% % 
% % 
% % 
% % % % % evolution of determinent
% % % figure(11);hold on
% % % plot(0:dt_store:Tmax,eigdet,'Color',color,'LineWidth',1)
% % % xlabel('Time (s)');ylabel('determinent of A','FontSize',14)
% % % set(gca,'Xtick',Tinit(21:20:length(Tinit)),'FontSize',14)
% % % set(gca,'XTickLabel',(20:20:length(Tinit)-1)*(Tstim+Tmemory+Tforget)/1000)
% % % xlim([Tinit(1)-500 Tmax-500])
% % 
% % % % weight change
% % figure(12);hold on
% % title('middle row of MEE')
% % plot(x,MEE0(nx/2,:),'g')
% % plot(x,MEEt(nx/2,:,1))
% % plot(x,MEEt(nx/2,:,end),'r')
% % legend('before perturbation','after perturbation','after learning')
% % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % xlim([-pi pi]);
% % hold off
% % % %
% % 
% % % Recurrent Input
% % RecInput_EE = MEE*SEE(:,end);
% % RecInput_EI = MEI*SEI(:,end);
% % figure(4);plot(x,RecInput_EE/max(RecInput_EE),'b','LineWidth',0.5);
% % hold on;
% % plot(x,RecInput_EI/max(RecInput_EE),'r','LineWidth',0.5);
% % plot(x(index),RecInput_EE(index)/max(RecInput_EE),'b','LineWidth',0.5,'Marker','o','MarkerSize',10);
% % plot(x(index),RecInput_EI(index)/max(RecInput_EE),'r','LineWidth',0.5,'Marker','o','MarkerSize',10);
% % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % xlim([-pi pi]);
% % ylim([0 1])
% % hold off
% % 
% % % % Spatial pattern of activity at every 1 s
% % figure(7);hold on
% % for l = 0:round(Tmemory/1000)
% %     plot(x,RE(:,(Tinit(end-1)+1000*l)/dt_store),'Color',color,'LineWidth',0.5)
% %     plot(x(index),RE(index,(Tinit(end-1)+1000*l)/dt_store),'Color',color,'LineWidth',0.5,'Marker','o','MarkerSize',10);
% % end
% % xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % xlim([-pi pi]);
% % 
% % set(gca,'Ytick',0:50:100,'FontSize',14)
% % ylim([0 120])
% % 
% % l = 3;
% % xdata = (round((Tinit(end-1)+Tstim+500)):dt_store:round(Tmax))';
% % ydata = RE(:,(Tinit(end-1)+1000*l)/dt_store)';
% % myfun = @(p,xdata) p(1)*exp(-(x-p(2)).^2/p(3))+p(4);
% % options=optimset('Display','off');
% % [p,resnorm] = lsqcurvefit(myfun,[1;0;1;10],xdata,ydata,[],[],options);
% % plot(x,p(1)*exp(-(x-p(2)).^2/p(3))+p(4),'k','LineWidth',0.5);
% % 
% % % kmax = 11;
% % % for i = 1:kmax;
% % % %     coeff_test(i) = dx/pi*sum(RE(:,Tinit+1000*l).*cos((i-1)*x)');
% % %     coeff_test(i) = dx/pi*sum(RE(:,(Tinit(end-1)+1000*l)/dt_store).*cos((i-1)*x)');
% % % end
% % % k = 2;
% % % coeff_test(k+1)/coeff_test(k+2)
% % %  exp(p(3)*((k+1)^2-k^2)/4)
% % % l = 2;
% % % xdata = x;
% % % ydata = RE(:,(Tinit+1000*l))';
% % % myfun = @(p,xdata) p(1)*exp(-(x-p(2)).^2/p(3))+p(4);
% % % options=optimset('NonlEqnAlgorithm','lm');
% % % [p,resnorm] = lsqcurvefit(myfun,[1;0;1;0],xdata,ydata,[],[],options);
% % % sqrt(p(3))
% % % hold off
% % % % 
% % % figure(8);hold on
% % % for l = 1:round(Tmemory/1000)
% % %     plot(x,RI(:,(Tinit+Tstim+1000*l)),'r','LineWidth',2);
% % % end
% % % xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% % % xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% % % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % % xlim([-pi pi]);
% % % set(gca,'Ytick',-10:10:50,'FontSize',14)
% % % ylim([0 50]);
% % % hold off
% % % % 
% % l = 3;
% % % T = 1000*l+Tinit+Tstim;
% % T = 1000*l+Tinit(end-1)+Tstim; T = T/dt_store;
% % figure(9);plot(x,(RE(:,T)-min(RE(:,T)))/(max(RE(:,T))-min(RE(:,T))),'b','LineWidth',0.5);
% % xdata = x;
% % ydata = RE(:,end)';
% % myfun = @(p,xdata) p(1)*exp(-(x-p(2)).^2/p(3))+p(4);
% % options=optimset('Display','off');
% % % options=optimset('NonlEqnAlgorithm','lm');
% % [p,resnorm] = lsqcurvefit(myfun,[1;0;1;0],xdata,ydata,[],[],options);
% % hold on;
% % plot(x,(p(1)*exp(-(x-p(2)).^2/p(3))+p(4)-min(RE(:,T)))/(max(RE(:,T))-min(RE(:,T))),'b','LineWidth',0.5);
% % 
% % plot(x,(RI(:,T)-min(RI(:,T)))/(max(RI(:,T))-min(RI(:,T))),'r','LineWidth',0.5)
% % 
% % ydata = RI(:,end)';
% % myfun = @(p,xdata) p(1)*exp(-(x-p(2)).^2/p(3))+p(4);
% % options=optimset('Display','off');
% % % options=optimset('NonlEqnAlgorithm','lm');
% % [p,resnorm] = lsqcurvefit(myfun,[1;0;1;20],xdata,ydata,[],[],options);
% % plot(x,(p(1)*exp(-(x-p(2)).^2/p(3))+p(4)-min(RI(:,T)))/(max(RI(:,T))-min(RI(:,T))),'r','LineWidth',0.5)
% % 
% % plot(x(index),(RE(index,T)-min(RE(index,T)))/(max(RE(index,T))-min(RE(index,T))),'b','LineWidth',0.5,'Marker','o','MarkerSize',10);
% % plot(x(index),(RI(index,T)-min(RI(index,T)))/(max(RI(index,T))-min(RI(index,T))),'r','LineWidth',0.5,'Marker','o','MarkerSize',10)
% % xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% % legend('Ex','In')
% % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % xlim([-pi pi]);
% % set(gca,'Ytick',0:0.5:1,'FontSize',14)
% % ylim([0 1]);
% % hold off
% % 
% % % % figure(10);hold on
% % % % for l = 1:round(Tmemory/1000)
% % % %     plot(x,Input_RE(:,(Tinit+Tstim+1000*l)/dt),'r','LineWidth',2);
% % % % end
% % % % xlabel('\theta','FontSize',14);ylabel('Firing Rate (Hz)','FontSize',14)
% % % % set(gca,'Xtick',pi*(-1:0.5:1),'FontSize',14)
% % % % set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
% % % % xlim([-pi pi]);
% % % % set(gca,'Ytick',-10:10:50,'FontSize',14)
% % % % ylim([0 50]);
% % % % hold off
% % 
%% save results
save([datapath,'/param.mat'],'-struct','param');
save([datapath,'/results.mat'],'t','TDelayOff','RE','RI','MEEt');