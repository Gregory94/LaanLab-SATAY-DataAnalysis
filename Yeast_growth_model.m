%% Growth of WT and bem1del /seperate experiments

% noc is number of colonies but measured as OD, tehrefore noc is here assumed to be directly related to the OD
% dt is doubling time
% WT is wild type
% del represents the bem1delete phenotype

runtime = 4800; %min
t = linspace(1,runtime,runtime);
start_OD = 0.01;
blank = 0.0085;
noc_WT = start_OD-blank; 
dt_WT = 78; 
noc_del = start_OD-blank; 
dt_del_mean = 229;
dt_del_min = 152;
dt_del_max = 301;
lag_WT = linspace(noc_WT,noc_WT,600);    %% possibility to assume different lagtimes for both strains
lag_del = linspace(noc_WT,noc_WT,600);

k_WT=log(2)/dt_WT;                      %% calculate growth constant from doubling time
k_del_mean=log(2)/dt_del_mean;
k_del_min=log(2)/dt_del_min;
k_del_max=log(2)/dt_del_max;

GC_WT = noc_WT*exp(k_WT*t);             %% caculate growth curve
GC_del_mean = noc_del*exp(k_del_mean*t);
GC_del_min = noc_del*exp(k_del_min*t);
GC_del_max = noc_del*exp(k_del_max*t);
GC_WT = cat(2,lag_WT,GC_WT);            %% add lag time to start
GC_del_mean = cat(2,lag_del,GC_del_mean);
GC_del_min = cat(2,lag_del,GC_del_min);
GC_del_max = cat(2,lag_del,GC_del_max);
GC_WT = GC_WT(1:runtime);               %% cut end to be the same length
GC_del_mean = GC_del_mean(1:runtime);
GC_del_min = GC_del_min(1:runtime);
GC_del_max = GC_del_max(1:runtime);

figure(1)
plot(t/60,GC_del_min,'LineWidth',2.0)
hold on
x2 = [t/60, fliplr(t/60)];
inBetween = [GC_del_min, fliplr(GC_del_max)];
fill(x2, inBetween, [0.9 0.95 0.95]);
hold on
plot(t/60,GC_del_max,'LineWidth',2.0)
hold on
plot(t/60,GC_del_mean,'LineWidth',2.0)
hold on
plot(t/60,GC_WT,'LineWidth',2.0)
set(gca,'FontSize',20)
%xlim([0,1000])
ylim([0,0.3])
legend({'bem1del-std','std','bem1del+std','bem1del','WT'},'Fontsize', 20)
xlabel('Time [h]','Fontsize', 30)
ylabel('OD600','Fontsize', 30)

%% Growth of WT and bem1del / single experiment

RWD = 0.00046; %Ratio WT to del, calculated as AID, bem3 and OsTIR1 length over the total yeast genome size

noc_WT1 = noc_WT*RWD;
noc_del1 = noc_del*(1-RWD);
lag_WT = linspace(noc_WT1,noc_WT1,600);    %% assuming a different lagtime for both strains
lag_del = linspace(noc_del1,noc_del1,1100);

GC_WT1 = noc_WT1*exp(k_WT*t);             %% caculate growth curve
GC_del1 = noc_del1*exp(k_del_mean*t);
GC_del1_min = noc_del1*exp(k_del_min*t);
GC_del1_max = noc_del1*exp(k_del_max*t);
GC_WT1 = cat(2,lag_WT,GC_WT1);            %% add lag time to start
GC_del1 = cat(2,lag_del,GC_del1);
GC_del1_min = cat(2,lag_del,GC_del1_min);
GC_del1_max = cat(2,lag_del,GC_del1_max);
GC_WT1 = GC_WT1(1:runtime);               %% cut end to be the same length
GC_del1 = GC_del1(1:runtime);
GC_del1_min = GC_del1_min(1:runtime);
GC_del1_max = GC_del1_max(1:runtime);

figure(1)
plot(t/60,GC_del1_min+blank,'LineWidth',2.0)
hold on
x2 = [t/60, fliplr(t/60)];
inBetween = [GC_del1_min+blank, fliplr(GC_del1_max+blank)];
fill(x2, inBetween, [0.9 0.95 0.95]);
hold on
plot(t/60,GC_del1_max+blank,'LineWidth',2.0)
hold on
plot(t/60,GC_del1+blank,'LineWidth',2.0)
hold on
plot(t/60,GC_WT1+blank,'LineWidth',2.0)
set(gca,'FontSize',20)
%xlim([0,1000])
ylim([0,0.3])
legend({'bem1del-std','std','bem1del+std','bem1del','WT'},'Fontsize', 20)
xlabel('time [h]','Fontsize', 30)
ylabel('OD600','Fontsize', 30)

rel_WT_min=GC_WT1./(GC_WT1+GC_del1_min);
rel_WT=GC_WT1./(GC_WT1+GC_del1);
rel_WT_max=GC_WT1./(GC_WT1+GC_del1_max);

figure(2)
plot(t/60, rel_WT_min,'LineWidth',2.0)
hold on
x2 = [t/60, fliplr(t/60)];
inBetween = [rel_WT_min, fliplr(rel_WT_max)];
fill(x2, inBetween, [0.9 0.95 0.95]);
hold on
plot(t/60, rel_WT_max,'LineWidth',2.0)
hold on
plot(t/60, rel_WT,'LineWidth',2.0)
set(gca,'FontSize',20)
legend({'[WT]/[WT]+[bem1del - std]','std','[WT]/[WT]+[bem1del + std]','[WT]/[WT]+[bem1del]'},'Fontsize', 30)
xlabel('time [h]','Fontsize', 30)
ylabel('[WT]/[WT]+[bem1del]','Fontsize', 30)