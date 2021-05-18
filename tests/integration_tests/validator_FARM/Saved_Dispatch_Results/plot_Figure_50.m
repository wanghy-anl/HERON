clear; tic
%% Specify the file name to read: No FARM
filename_1 = 'ERCOT_2_no_FARM_3Units';
fid = fopen(filename_1);
% Get one line from the file
tline = fgetl(fid);
% when this line is not empty (the end of file)
time_1=[];
BOP_vyminmax_1=[];
SES_vyminmax_1=[];
TES_vyminmax_1=[];
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
            BOP_vyminmax_1=[BOP_vyminmax_1;data];
            time_1 = [time_1;t_temp];
        elseif startsWith(tline, "SES ,")
            SES_vyminmax_1=[SES_vyminmax_1;data];
        elseif startsWith(tline, "TES ,")
            TES_vyminmax_1=[TES_vyminmax_1;data];
        end
        
    end
%     disp(tline)
    tline = fgetl(fid);
end
fclose(fid);

power_provided_1=BOP_vyminmax_1(:,1)+SES_vyminmax_1(:,1)+TES_vyminmax_1(:,1);
time_hour_1=[];power_array_hour_1=[];
for i=1:numel(time_1)
    if mod(time_1(i),3600)==1800
        time_hour_1 = [time_hour_1; time_1(i)];
        power_array_hour_1=[power_array_hour_1;BOP_vyminmax_1(i,1) SES_vyminmax_1(i,1) TES_vyminmax_1(i,1)];
    end
end
% convert output power to MW, convert output pressure to Bar
for i=2:4
    BOP_vyminmax_1(:,i)=BOP_vyminmax_1(:,i)*1e-6;
    SES_vyminmax_1(:,i)=SES_vyminmax_1(:,i)*1e-6;
end
for i=5:7
    BOP_vyminmax_1(:,i)=BOP_vyminmax_1(:,i)*1e-5;
end
time_hour_1 = time_hour_1/3600;
time_1 = time_1/3600;

%% Specify the file name to read: with FARM
filename_2 = 'ERCOT_1_FARM_3Units';
fid = fopen(filename_2);
% Get one line from the file
tline = fgetl(fid);
% when this line is not empty (the end of file)
time_2=[];
BOP_vyminmax_2=[];
SES_vyminmax_2=[];
TES_vyminmax_2=[];
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
            BOP_vyminmax_2=[BOP_vyminmax_2;data];
            time_2 = [time_2;t_temp];
        elseif startsWith(tline, "SES ,")
            SES_vyminmax_2=[SES_vyminmax_2;data];
        elseif startsWith(tline, "TES ,")
            TES_vyminmax_2=[TES_vyminmax_2;data];
        end
        
    end
%     disp(tline)
    tline = fgetl(fid);
end
fclose(fid);

power_provided_2=BOP_vyminmax_2(:,1)+SES_vyminmax_2(:,1)+TES_vyminmax_2(:,1);
time_hour_2=[];power_array_hour_2=[];
for i=1:numel(time_2)
    if mod(time_2(i),3600)==1800
        time_hour_2 = [time_hour_2; time_2(i)];
        power_array_hour_2=[power_array_hour_2;BOP_vyminmax_2(i,1) SES_vyminmax_2(i,1) TES_vyminmax_2(i,1)];
    end
end
% convert output power to MW, convert output pressure to Bar
for i=2:4
    BOP_vyminmax_2(:,i)=BOP_vyminmax_2(:,i)*1e-6;
    SES_vyminmax_2(:,i)=SES_vyminmax_2(:,i)*1e-6;
end
for i=5:7
    BOP_vyminmax_2(:,i)=BOP_vyminmax_2(:,i)*1e-5;
end
time_hour_2 = time_hour_2/3600;
time_2 = time_2/3600;


%% 4. Plot the implicit constraints in 2D 
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.05], [0.1 0.05], [0.08 0.05]);

figure(50)
set(gcf,'Position',[100 100 1400 600])

for unit_idx=1:3
    % Plot y1(x) with y2 (y) and their min/max
    subplot(1,3,unit_idx)
    if unit_idx==1 % BOP
        x_lb = min(BOP_vyminmax_1(:,3)); x_ub = max(BOP_vyminmax_1(:,4));
        y_lb = min(BOP_vyminmax_1(:,6)); y_ub = max(BOP_vyminmax_1(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(BOP_vyminmax_1(:,2),BOP_vyminmax_1(:,5),'-o','MarkerSize',8)
        plot(BOP_vyminmax_2(:,2),BOP_vyminmax_2(:,5),'-^','LineWidth',1)       
        hold off        
        title("BOP Implicit Constraints"); 
        ylabel('y2, Turbine Pressure (bar)'); 
        xlabel('y1, Output Power (MW)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(a)','FontSize',32)
        legend('Without FARM', 'With FARM')
    elseif unit_idx==2
        x_lb = min(SES_vyminmax_1(:,3)); x_ub = max(SES_vyminmax_1(:,4));
        y_lb = min(SES_vyminmax_1(:,6)); y_ub = max(SES_vyminmax_1(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(SES_vyminmax_1(:,2),SES_vyminmax_1(:,5),'-o','MarkerSize',8)
        plot(SES_vyminmax_2(:,2),SES_vyminmax_2(:,5),'-^','LineWidth',1)
        hold off
        title("SES Implicit Constraints"); 
        ylabel('y2, Firing Temperature (K)'); 
        xlabel('y1, Output Power (MW)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(b)','FontSize',32)
        legend('Without FARM', 'With FARM')
    elseif unit_idx==3
        x_lb = min(TES_vyminmax_1(:,3)); x_ub = max(TES_vyminmax_1(:,4));
        y_lb = min(TES_vyminmax_1(:,6)); y_ub = max(TES_vyminmax_1(:,7));
        rectangle('Position',[x_lb y_lb x_ub-x_lb y_ub-y_lb],'LineStyle','--','LineWidth',3)
        xlim([x_lb-(x_ub-x_lb)*0.2 x_ub+(x_ub-x_lb)*0.2])
        ylim([y_lb-(y_ub-y_lb)*0.2 y_ub+(y_ub-y_lb)*0.2])
        hold on
        plot(TES_vyminmax_1(:,2),TES_vyminmax_1(:,5),'-o','MarkerSize',8)
        plot(TES_vyminmax_2(:,2),TES_vyminmax_2(:,5),'-^','LineWidth',1)
        hold off
        title("TES Implicit Constraints"); 
        ylabel('y2, Cold Tank Level (m)') 
        xlabel('y1, Hot Tank Level (m)')
%         text(x_lb-(x_ub-x_lb)*0.1, y_ub+(y_ub-y_lb)*0.1, '(c)','FontSize',32)
        legend('Without FARM', 'With FARM')
    end
end
print('Figure_50.png','-dpng','-r300')



%%
toc
