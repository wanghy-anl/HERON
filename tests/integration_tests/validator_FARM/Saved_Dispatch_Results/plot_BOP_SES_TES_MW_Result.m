clear; tic

% Specify the file name to read
filename = 'ERCOT_1_FARM_3Units'; title_appendix = ', with FARM';
% filename = 'ERCOT_2_no_FARM_3Units'; title_appendix = ', without FARM';
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
%% 1. Plot the power dispatch stack
x_label_min=0;
x_label_max=168;
x_tick_interval=24;

figure(10)
set(gcf,'Position',[100 50 600 500])
% Plot the stacked bar of power components
bar(time_hour, power_array_hour, 'stacked');hold on
% Plot the total power provided
plot(time, power_provided,'LineWidth',3);hold off
xlabel('Time (Hour)');ylabel('Power (MW)'); 
xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
ylim([-400 1600])
% legend('BOP Output Power','SES Output Power','TES Discharging(+)/Charging(-)','Market Demand','Location','best')
legend('BOP','SES','TES Discharging(+)/Charging(-)','Market Demand','Location','best')
% title(strcat('Contribution of each Power Source',title_appendix))

print(strrep(strrep(strcat('Figure_10',title_appendix,'.png'),',','_'),' ','_'),'-dpng','-r300')

%% 2. Plot the explicit and implicit constraints v.s. time
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.1 0.05], [0.08 0.05]);

figure(20)
set(gcf,'Position',[100 50 1400 1050])

for unit_idx=1:3
    % Plot the output power for all 3 units
    subplot(3,3,(unit_idx-1)*3+1)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,1))
%         y_lb = min(BOP_vyminmax(:,1)); y_ub = max(BOP_vyminmax(:,1));
        y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
        title(strcat('BOP Dispatched Power', title_appendix))
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(a)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,1))
%         y_lb = min(SES_vyminmax(:,1)); y_ub = max(SES_vyminmax(:,1));
        y_lb = -5; y_ub = max(SES_vyminmax(:,4));
        title(strcat('SES Dispatched Power', title_appendix))
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(d)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,1))
        y_lb = min(TES_vyminmax(:,1)); y_ub = max(TES_vyminmax(:,1));
        title({strcat('TES Dispatched Power', title_appendix),'Discharging(+)/Charging(-)'})
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(g)','FontSize',32)
    end
    xlabel('Time (Hour)');ylabel('Power (MW)'); 
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    % TODO: Change the scale to the y1 level for BOP and SES
%     ylim([y_lb-1 y_ub+1])
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
%     ytickformat('%.2f')
    
    % Plot y1 and its min/max
    subplot(3,3,(unit_idx-1)*3+2)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,4),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,2),'-k'); % y1
        plot(time,BOP_vyminmax(:,3),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
        title("BOP Constraint 1: Output Power"); ylabel('Power (MW)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(b)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,2),'-k');
        plot(time,SES_vyminmax(:,3),'--b','LineWidth',3); hold off
%         y_lb = min(SES_vyminmax(:,3)); y_ub = max(SES_vyminmax(:,4));
        y_lb = -5; y_ub = max(SES_vyminmax(:,4));
        title("SES Constraint 1: Output Power"); ylabel('Power (MW)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(e)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,2),'-k');
        plot(time,TES_vyminmax(:,3),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,3)); y_ub = max(TES_vyminmax(:,4));
        title("TES Constraint 1: Hot Tank Level"); ylabel('Level (m)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(h)','FontSize',32)
    end
    xlabel('Time (Hour)');
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #1','Lower Bound','Location','southwest')
    
    % Plot y2 and its min/max
    subplot(3,3,(unit_idx-1)*3+3)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,7),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,5),'-k'); % y1
        plot(time,BOP_vyminmax(:,6),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,6)); y_ub = max(BOP_vyminmax(:,7));
        title("BOP Constraint 2: Turbine Pressure"); ylabel('Pressure (Bar)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(c)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,5),'-k');
        plot(time,SES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(SES_vyminmax(:,6)); y_ub = max(SES_vyminmax(:,7));
        title("SES Constraint 2: Firing Temperature"); ylabel('Temperature (K)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(f)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,5),'-k');
        plot(time,TES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,6)); y_ub = max(TES_vyminmax(:,7));
        title("TES Constraint 2: Cold Tank Level"); ylabel('Level (m)'); 
%         text(12, y_ub-(y_ub-y_lb)*0.13, '(i)','FontSize',32)
    end
    xlabel('Time (Hour)');
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #2','Lower Bound','Location','southwest')
        
    
end


print(strrep(strrep(strcat('Figure_20',title_appendix,'.png'),',','_'),' ','_'),'-dpng','-r300')

%% 3. Plot everything in one figure
figure(30)
set(gcf,'Position',[100 50 2240 1260])
subplot(3,4,[1 5 9]) 
% Plot the stacked bar of power components
bar(time_hour, power_array_hour, 'stacked');hold on
% Plot the total power provided
plot(time, power_provided,'LineWidth',3);hold off
xlabel('Time (Hour)');ylabel('Power (MW)'); 
xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
legend('BOP Output Power','SES Output Power','TES Discharging(+)/Charging(-)','Market Demand','Location','best')
title(strcat('Contribution of each Power Source',title_appendix))

for unit_idx=1:3
    % Plot the output power for all 3 units
    subplot(3,4,(unit_idx-1)*4+2)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,1))
%         y_lb = min(BOP_vyminmax(:,1)); y_ub = max(BOP_vyminmax(:,1));
        y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
        title(strcat('BOP Dispatched Power', title_appendix))
        text(12, y_ub-(y_ub-y_lb)*0.13, '(a)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,1))
%         y_lb = min(SES_vyminmax(:,1)); y_ub = max(SES_vyminmax(:,1));
        y_lb = -5; y_ub = max(SES_vyminmax(:,4));
        title(strcat('SES Dispatched Power', title_appendix))
        text(12, y_ub-(y_ub-y_lb)*0.13, '(d)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,1))
        y_lb = min(TES_vyminmax(:,1)); y_ub = max(TES_vyminmax(:,1));
        title({strcat('TES Dispatched Power', title_appendix),'Discharging(+)/Charging(-)'})
        text(12, y_ub-(y_ub-y_lb)*0.13, '(g)','FontSize',32)
    end
    xlabel('Time (Hour)');ylabel('Power (MW)'); 
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    % TODO: Change the scale to the y1 level for BOP and SES
%     ylim([y_lb-1 y_ub+1])
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
%     ytickformat('%.2f')
    
    % Plot y1 and its min/max
    subplot(3,4,(unit_idx-1)*4+3)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,4),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,2),'-k'); % y1
        plot(time,BOP_vyminmax(:,3),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,3)); y_ub = max(BOP_vyminmax(:,4));
        title("BOP Constraint 1: Output Power"); ylabel('Power (MW)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(b)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,2),'-k');
        plot(time,SES_vyminmax(:,3),'--b','LineWidth',3); hold off
%         y_lb = min(SES_vyminmax(:,3)); y_ub = max(SES_vyminmax(:,4));
        y_lb = -5; y_ub = max(SES_vyminmax(:,4));
        title("SES Constraint 1: Output Power"); ylabel('Power (MW)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(e)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,4),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,2),'-k');
        plot(time,TES_vyminmax(:,3),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,3)); y_ub = max(TES_vyminmax(:,4));
        title("TES Constraint 1: Hot Tank Level"); ylabel('Level (m)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(h)','FontSize',32)
    end
    xlabel('Time (Hour)');
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #1','Lower Bound','Location','southwest')

    % Plot y2 and its min/max
    subplot(3,4,(unit_idx-1)*4+4)
    if unit_idx==1
        plot(time,BOP_vyminmax(:,7),'--r','LineWidth',3); hold on % y1max
        plot(time,BOP_vyminmax(:,5),'-k'); % y1
        plot(time,BOP_vyminmax(:,6),'--b','LineWidth',3); hold off %y1min
        y_lb = min(BOP_vyminmax(:,6)); y_ub = max(BOP_vyminmax(:,7));
        title("BOP Constraint 2: Turbine Pressure"); ylabel('Pressure (Bar)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(c)','FontSize',32)
    elseif unit_idx==2
        plot(time,SES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,SES_vyminmax(:,5),'-k');
        plot(time,SES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(SES_vyminmax(:,6)); y_ub = max(SES_vyminmax(:,7));
        title("SES Constraint 2: Firing Temperature"); ylabel('Temperature (K)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(f)','FontSize',32)
    elseif unit_idx==3
        plot(time,TES_vyminmax(:,7),'--r','LineWidth',3); hold on
        plot(time,TES_vyminmax(:,5),'-k');
        plot(time,TES_vyminmax(:,6),'--b','LineWidth',3); hold off
        y_lb = min(TES_vyminmax(:,6)); y_ub = max(TES_vyminmax(:,7));
        title("TES Constraint 2: Cold Tank Level"); ylabel('Level (m)'); 
        text(12, y_ub-(y_ub-y_lb)*0.13, '(i)','FontSize',32)
    end
    xlabel('Time (Hour)');
    xlim([x_label_min x_label_max]);xticks(x_label_min:x_tick_interval:x_label_max)
    ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
    legend('Upper Bound','Output #2','Lower Bound','Location','southwest')
        
    
end
print(strrep(strrep(strcat('Figure_30',title_appendix,'.png'),',','_'),' ','_'),'-dpng','-r300')

%% 4. Plot the implicit constraints in 2D 
figure(40)
set(gcf,'Position',[100 100 1400 600])

for unit_idx=1:3
    % Plot y1(x) with y2 (y) and their min/max
    subplot(1,3,unit_idx)
    if unit_idx==1 % BOP
        x_lb = min(BOP_vyminmax(:,3)); x_ub = max(BOP_vyminmax(:,4));
        y_lb = min(BOP_vyminmax(:,6)); y_ub = max(BOP_vyminmax(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(BOP_vyminmax(:,2),BOP_vyminmax(:,5),'-o')
        hold off        
        title(strcat('BOP Implicit Constraints',title_appendix))
        ylabel('y2, Turbine Pressure (Bar)'); 
        xlabel('y1, Output Power (MW)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(a)','FontSize',32)
%         legend('Turbine Pressure v.s. Output Power')
    elseif unit_idx==2
        x_lb = min(SES_vyminmax(:,3)); x_ub = max(SES_vyminmax(:,4));
        y_lb = min(SES_vyminmax(:,6)); y_ub = max(SES_vyminmax(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(SES_vyminmax(:,2),SES_vyminmax(:,5),'-o')
        hold off
        title(strcat('SES Implicit Constraints',title_appendix))
        ylabel('y2, Firing Temperature (K)'); 
        xlabel('y1, Output Power (MW)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(b)','FontSize',32)
%         legend('Firing Temperature v.s. Output Power')
    elseif unit_idx==3
        x_lb = min(TES_vyminmax(:,3)); x_ub = max(TES_vyminmax(:,4));
        y_lb = min(TES_vyminmax(:,6)); y_ub = max(TES_vyminmax(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(TES_vyminmax(:,2),TES_vyminmax(:,5),'-o')
        hold off
        title(strcat('TES Implicit Constraints',title_appendix))
        ylabel('y2, Cold Tank Level (m)') 
        xlabel('y1, Hot Tank Level (m)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(c)','FontSize',32)
%         legend('Cold Tank Level v.s. Hot Tank Level')
    end
end
print(strrep(strrep(strcat('Figure_40',title_appendix,'.png'),',','_'),' ','_'),'-dpng','-r300')



%%
toc