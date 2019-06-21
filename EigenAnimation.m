param = load('param.mat');
x = param.x;
load('results.mat','MEEt');
figure('Position',[50 100 1250 500])
nn=5;c = winter(nn);%c=fliplr(c);
[V,D] = sorteigs((MEEt(:,:,1)-eye(param.N)) - param.MEI/ (param.MII + eye(param.N)) * param.MIE);
subplot(1,2,1)
hold on
line1 = plot(real(diag(D)));
points1 = scatter(1:nn,real(diag(D(1:nn,1:nn))),80,c,'filled');
ylim([-3.5,.5])
ylabel('real(\lambda)')
t1=title(sprintf('%d ms',param.TDelayOff(1)));
hold off
subplot(1,2,2)
hold on
for k=1:nn
    line2(k) = plot(x,real(V(:,k)),'Color',c(k,:));
end
legend(num2str((1:nn)','No.%d eigenvector'),'Location','southeast')
xlabel('\theta');
set(gca,'Xtick',pi*(-1:0.5:1))
set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
xlim([-pi pi]);
ylabel('real(v)')
ylim([-.2,.2])
title('eigenvectors,sorted by largest fourier coff')
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(:,:,1,1) = rgb2ind(f.cdata,map,'nodither');
for i=2:size(MEEt,3)
    [V,D] = sorteigs((MEEt(:,:,i)-eye(param.N)) - param.MEI / (param.MII + eye(param.N)) * param.MIE);
%     line1.YData = real(diag(D));
    set(line1,'YData',real(diag(D)));
%     points1.YData = real(diag(D(1:nn,1:nn)));
    set(points1,'Ydata',real(diag(D(1:nn,1:nn))));
    ylim([-3.5,.5])
    ylabel('real(\lambda)')
%     t1.String=sprintf('eigenvalues %d s',i*dt_store/1000);
    set(t1,'String',sprintf('eigenvalues %d ms',param.TDelayOff(i)));
    subplot(1,2,2)
    for k=1:nn
%         line2(k).YData=abs(V(:,k));
        set(line2(k),'YData',real(V(:,k)));
    end

    ylabel('real(v)')
    ylim([-.2,.2])
    f = getframe(gcf);
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'eigenvalues_vectors.gif','DelayTime',0,'LoopCount',inf)
