% Plot the ERCOT comparison

% Specify the file name
xls_filename = 'netLoad_Comparison.xlsx';
T = readtable(xls_filename);
T = T{:,:};

figure(1)
set(gcf,'Position',[500 400 1000 450])
plot(T(:,1),-1*T(:,2),'-o','Color','#D95319','LineWidth',2,'MarkerEdgeColor','#D95319','MarkerSize',6)
hold on
plot(T(:,1),-1*T(:,3),'-^','Color','#0072BD','LineWidth',2,'MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD','MarkerSize',5)
hold off
legend('ERCOT Training Data', 'Fourier/ARMA Modeled Data')
xlabel('Time (Hours)')
ylabel('Power (MW)')

x_label_min=0;
x_label_max=168;
x_tick_interval=24;
xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
ylim([1080 1300])
grid on
set(gca,'FontSize',12)
%%
print('ERCOT_Comparison.png','-dpng','-r300')
