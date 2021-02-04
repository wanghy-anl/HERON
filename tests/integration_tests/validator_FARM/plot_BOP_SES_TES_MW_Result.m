clear; tic
% Specify the file name to read
filename = 'Sweep_Runs_o\sweep\1\out~inner';
fid = fopen(filename);
% Get one line from the file
tline = fgetl(fid);
% when this line is not empty (the end of file)
time=[];
BOP_vyminmax=[];
SES_vyminmax=[];
TES_vyminmax=[];
while ischar(tline)
    if any([startsWith(tline, "BOP ,") startsWith(tline, "SES ,") startsWith(tline, "TES ,")]) 
        data=nan(1,7);
        c = strsplit(tline,',');
        for i=1:numel(c)
            if strcmp(c{i},'t')
                t_temp = str2double(c{i+1});
            elseif strcmp(c{i},'vp')
                data(1)=str2double(c{i+1});
            elseif strcmp(c{i},'y1')
                data(2)=str2double(c{i+1});
            elseif strcmp(c{i},'y1min')
                data(3)=str2double(c{i+1});
            elseif strcmp(c{i},'y1max')
                data(4)=str2double(c{i+1});
            elseif strcmp(c{i},'y2')
                data(5)=str2double(c{i+1});
            elseif strcmp(c{i},'y2min')
                data(6)=str2double(c{i+1});
            elseif strcmp(c{i},'y2max')
                data(7)=str2double(c{i+1});
            end
        end
        
        if startsWith(tline, "BOP ,")
            BOP_vyminmax=[BOP_vyminmax;data];
            time = [time;t_temp];
        elseif startsWith(tline, "SES ,")
            SES_vyminmax=[SES_vyminmax;data];
        elseif startsWith(tline, "TES ,")
            TES_vyminmax=[TES_vyminmax;data];
        end
        
    end
%     disp(tline)
    tline = fgetl(fid);
end
fclose(fid);

power_provided=BOP_vyminmax(:,1)+SES_vyminmax(:,1)+TES_vyminmax(:,1);
time_hour=[];power_array_hour=[];
for i=1:numel(time)
    if mod(time(i),3600)==1800
        time_hour = [time_hour; time(i)];
        power_array_hour=[power_array_hour;BOP_vyminmax(i,1) SES_vyminmax(i,1) TES_vyminmax(i,1)];
    end
end
% convert output power to MW, convert output pressure to Bar
for i=2:4
    BOP_vyminmax(:,i)=BOP_vyminmax(:,i)*1e-6;
    SES_vyminmax(:,i)=SES_vyminmax(:,i)*1e-6;
end
for i=5:7
    BOP_vyminmax(:,i)=BOP_vyminmax(:,i)*1e-5;
end
time_hour = time_hour/3600;
time = time/3600;
%%
% figure(1)
% set(gcf,'Position',[100 50 600 500])
% % Plot the stacked bar of power components
% bar(time_hour, power_array_hour, 'stacked');hold on
% % Plot the total power provided
% plot(time, power_provided,'LineWidth',3);hold off
% xlabel('Time (Hour)');ylabel('Power (MW)'); 
% xlim([0 24]);xticks(0:4:24)
% legend('BOP Output Power','SES Output Power','TES Discharging(+)/Charging(-)','Market Demand','Location','best')
% title('Contribution of each Power Source')
%%
% figure(2)
% set(gcf,'Position',[100 50 1400 1050])
% 
% for unit_idx=1:3
%     % Plot the output power for all 3 units
%     subplot(3,3,(unit_idx-1)*3+1)
%     if unit_idx==1
%         plot(time,BOP_vyminmax(:,1))
%         y_lb = min(BOP_vyminmax(:,1)); y_ub = max(BOP_vyminmax(:,1));
%         title('BOP Dispatched Power')
%     elseif unit_idx==2
%         plot(time,SES_vyminmax(:,1))
%         y_lb = min(SES_vyminmax(:,1)); y_ub = max(SES_vyminmax(:,1));
%         title('SES Dispatched Power')
%     elseif unit_idx==3
%         plot(time,TES_vyminmax(:,1))
%         y_lb = min(TES_vyminmax(:,1)); y_ub = max(TES_vyminmax(:,1));
%         title({'TES Dispatched Power','Discharging(-)/Charging(+)'})
%     end
%     xlabel('Time (Hour)');ylabel('Power (MW)'); 
%     xlim([0 24]); ylim([y_lb-1 y_ub+1])
%     xticks(0:4:24)
%     ytickformat('%.2f')
%     
%     % Plot y1 and its min/max
%     subplot(3,3,(unit_idx-1)*3+2)
%     if unit_idx==1
%         plot(time,BOP_vyminmax(:,4),'--r','LineWidth',3); hold on % y1max
%         plot(time,BOP_vyminmax(:,2),'-k'); % y1
%         plot(time,BOP_vyminmax(:,3),'--b','LineWidth',3); hold off %y1min
%         y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
%         title("BOP Constraint 1: Output Power"); ylabel('Power (MW)'); 
%     elseif unit_idx==2
%         plot(time,SES_vyminmax(:,4),'--r','LineWidth',3); hold on
%         plot(time,SES_vyminmax(:,2),'-k');
%         plot(time,SES_vyminmax(:,3),'--b','LineWidth',3); hold off
%         y_lb = min(SES_vyminmax(:,3)); y_ub = max(SES_vyminmax(:,4));
%         title("SES Constraint 1: Output Power"); ylabel('Power (MW)'); 
%     elseif unit_idx==3
%         plot(time,TES_vyminmax(:,4),'--r','LineWidth',3); hold on
%         plot(time,TES_vyminmax(:,2),'-k');
%         plot(time,TES_vyminmax(:,3),'--b','LineWidth',3); hold off
%         y_lb = min(TES_vyminmax(:,3)); y_ub = max(TES_vyminmax(:,4));
%         title("TES Constraint 1: Hot Tank Level"); ylabel('Level (m)'); 
%     end
%     xlabel('Time (Hour)');xlim([0 24]); xticks(0:4:24)
%     ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
%     legend('Upper Bound','Output #1','Lower Bound','Location','southwest')
%     
%     % Plot y2 and its min/max
%     subplot(3,3,(unit_idx-1)*3+3)
%     if unit_idx==1
%         plot(time,BOP_vyminmax(:,7),'--r','LineWidth',3); hold on % y1max
%         plot(time,BOP_vyminmax(:,5),'-k'); % y1
%         plot(time,BOP_vyminmax(:,6),'--b','LineWidth',3); hold off %y1min
%         y_lb = min(BOP_vyminmax(:,6)); y_ub = max(BOP_vyminmax(:,7));
%         title("BOP Constraint 2: Turbine Pressure"); ylabel('Pressure (Bar)'); 
%     elseif unit_idx==2
%         plot(time,SES_vyminmax(:,7),'--r','LineWidth',3); hold on
%         plot(time,SES_vyminmax(:,5),'-k');
%         plot(time,SES_vyminmax(:,6),'--b','LineWidth',3); hold off
%         y_lb = min(SES_vyminmax(:,6)); y_ub = max(SES_vyminmax(:,7));
%         title("SES Constraint 2: Firing Temperature"); ylabel('Temperature (K)'); 
%     elseif unit_idx==3
%         plot(time,TES_vyminmax(:,7),'--r','LineWidth',3); hold on
%         plot(time,TES_vyminmax(:,5),'-k');
%         plot(time,TES_vyminmax(:,6),'--b','LineWidth',3); hold off
%         y_lb = min(TES_vyminmax(:,6)); y_ub = max(TES_vyminmax(:,7));
%         title("TES Constraint 2: Cold Tank Level"); ylabel('Level (m)'); 
%     end
%     xlabel('Time (Hour)');xlim([0 24]); xticks(0:4:24)
%     ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
%     legend('Upper Bound','Output #2','Lower Bound','Location','southwest')
%         
%     
% end
%% Plot everything in one figure
figure(3)
set(gcf,'Position',[100 50 2240 1260])
subplot(3,4,[1 5 9]) 
% Plot the stacked bar of power components
bar(time_hour, power_array_hour, 'stacked');hold on
% Plot the total power provided
plot(time, power_provided,'LineWidth',3);hold off
xlabel('Time (Hour)');ylabel('Power (MW)'); 
xlim([0 24]);xticks(0:4:24)
legend('BOP Output Power','SES Output Power','TES Discharging(+)/Charging(-)','Market Demand','Location','best')
title('Contribution of each Power Source')

for unit_idx=1:3
    % Plot the output power for all 3 units
    subplot(3,4,(unit_idx-1)*4+2)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,1))
        y_lb = min(BOP_vyminmax(:,1)); y_ub = max(BOP_vyminmax(:,1));
        title('BOP Dispatched Power')
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,1))
        y_lb = min(SES_vyminmax(:,1)); y_ub = max(SES_vyminmax(:,1));
        title('SES Dispatched Power')
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,1))
        y_lb = min(TES_vyminmax(:,1)); y_ub = max(TES_vyminmax(:,1));
        title({'TES Dispatched Power','Discharging(-)/Charging(+)'})
    end
    xlabel('Time (Hour)');ylabel('Power (MW)'); 
    xlim([0 24]); ylim([y_lb-1 y_ub+1])
    xticks(0:4:24)
    ytickformat('%.2f')
    
    % Plot y1 and its min/max
    subplot(3,4,(unit_idx-1)*4+3)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,4),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,2),'-k'); % y1
        plot(time,BOP_vyminmax(:,3),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
        title("BOP Constraint 1: Output Power"); ylabel('Power (MW)'); 
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,2),'-k');
        plot(time,SES_vyminmax(:,3),'--b','LineWidth',3); hold off
        y_lb = min(SES_vyminmax(:,3)); y_ub = max(SES_vyminmax(:,4));
        title("SES Constraint 1: Output Power"); ylabel('Power (MW)'); 
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,2),'-k');
        plot(time,TES_vyminmax(:,3),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,3)); y_ub = max(TES_vyminmax(:,4));
        title("TES Constraint 1: Hot Tank Level"); ylabel('Level (m)'); 
    end
    xlabel('Time (Hour)');xlim([0 24]); xticks(0:4:24)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #1','Lower Bound','Location','best')
    % Plot y2 and its min/max
    subplot(3,4,(unit_idx-1)*4+4)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,7),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,5),'-k'); % y1
        plot(time,BOP_vyminmax(:,6),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,6)); y_ub = max(BOP_vyminmax(:,7));
        title("BOP Constraint 2: Turbine Pressure"); ylabel('Pressure (Bar)'); 
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,5),'-k');
        plot(time,SES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(SES_vyminmax(:,6)); y_ub = max(SES_vyminmax(:,7));
        title("SES Constraint 2: Firing Temperature"); ylabel('Temperature (K)'); 
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,5),'-k');
        plot(time,TES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,6)); y_ub = max(TES_vyminmax(:,7));
        title("TES Constraint 2: Cold Tank Level"); ylabel('Level (m)'); 
    end
    xlabel('Time (Hour)');xlim([0 24]); xticks(0:4:24)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #2','Lower Bound','Location','best')
        
    
end


%%
toc
