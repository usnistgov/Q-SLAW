
%% basic plots

frame_num = 5;

f1 = figure;
f1.Position = [100,100,1500,1000];
subplot(2,1,1)
u_curr = u{frame_num}{1}(1:end-0,1:end-0);
contourf(u_curr,20,'linestyle','none')
c= colorbar;
c.Label.String = 'Displacement (px)';
axis image
colormap(cmrMap)
xlabel('x_1 (px)')
ylabel('x_2 (px)')
title('u_1')
set(gca,'fontsize',24)

subplot(2,1,2)
u_curr = u{frame_num}{2}(1:end-0,1:end-0);
contourf(u_curr,20,'linestyle','none')
c = colorbar;
c.Label.String = 'Displacement (px)';
axis image
colormap(cmrMap)
xlabel('x_1 (px)')
ylabel('x_2 (px)')
title('u_2')
set(gca,'fontsize',24)

%saveas(gcf,'./VN01_001_003_QST1_002_displacement_frame06.png')

 
%% plot with median removed

figure
subplot(1,2,1)
u_curr = u{end}{1}(9:end-8,9:end-8);
contourf(u_curr - median(u_curr),'all','omitnan'),'linestyle','none')
c= colorbar;
c.Label.String = 'Displacement (px)';
axis image
colormap(cmrMap)
xlabel('x_1 (px)')
ylabel('x_2 (px)')
title('Displacement in the x_1 direction, u_1')
set(gca,'fontsize',24)

subplot(1,2,2)
u_curr = u{end}{2}(9:end-8,9:end-8);
contourf(u_curr - median(u_curr,'all','omitnan'),'linestyle','none')
c= colorbar;
c.Label.String = 'Displacement (px)';
axis image
colormap(cmrMap)
xlabel('x_1 (px)')
ylabel('x_2 (px)')
title('Displacement in the x_2 direction, u_2')
set(gca,'fontsize',24)



