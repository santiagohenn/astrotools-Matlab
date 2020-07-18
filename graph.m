figure
for k=1:nSM
subplot(1,nSM,k)
plot(1:1:schedule,dataSet(1:schedule,6,k),'--r.')
hold on
plot(1:1:schedule,dataSet(1:schedule,7,k),'--b.')
hold on
plot(1:1:schedule,dataSet(1:schedule,9,k),'--g.')
hold on
ylabel('Grados')
xlabel('Iteration')
grid on
legend ('Inclinacion','RAAN','Anomalía v.')
end

%%
figure
for k=1:nRelay
subplot(1,nRelay,k)
k=k+nMS;
plot(1:1:schedule,dataSet(:,6,k),'--r.')
hold on
plot(1:1:schedule,dataSet(:,7,k),'--b.')
hold on
plot(1:1:schedule,dataSet(:,8,k),'--g.')
hold on
ylabel('Grados')
xlabel('Iteration')
grid on
legend ('Inclinacion','RAAN','Arg. del Perigeo')
end

%% Plot profit for GA

figure('Position',[70 70 1000 400])
plot(1:1:generations,maxFitness,'Color',[0 0.547 0.841],'LineWidth',2); hold on;
plot(1:1:generations,totalFitness/population,'Color',[0.4 0.4 0.4],'LineWidth',2); hold on;
ylabel('Métrica [%]');
xlabel('Generación');
legend('Mejor Metrica hallada','Metrica promedio');
legend('Location','northeastoutside');
grid on;
set(gca,'GridAlpha',0.25)
annotation('textbox',[.725 .66 .1 .1],'String',strcat('Máximo valor hallado:',{' '},num2str(max(maxFitness))));
set(gca,'FontSize',12)
set(gca,'Position',[0.10038 0.158 0.59843 0.797])

%% Plot profit function for SA & BA & HC

figure('Position',[70 70 1000 400]) %[0 0.547 0.841]
plot(1:1:schedule,dataSet(1:schedule,1),'Color',[0.4 0.4 0.4],'LineWidth',2); hold on;
plot(1:1:schedule,dataSet(1:schedule,2),'Color',[0 0.547 0.841],'LineWidth',2); hold on;
plot(1:1:schedule,dataSet(1:schedule,3),'r','LineWidth',2); hold on;
ylabel('Métrica [%]');
xlabel('Iteración');
legend('Solución evaluada','Solución aceptada','Mejor solución hallada');
legend('Location','northeastoutside')
annotation('textbox',[.725 .66 .1 .1],'String',strcat('Máximo valor hallado:',{' '},num2str(max(dataSet(1:schedule,3)))));
set(gca,'FontSize',12)
set(gca,'GridAlpha',0.25)
set(gca,'Position',[0.10038 0.158 0.59843 0.797])
grid on
hold on

%% Plot profit function for SA & BA & HC

figure('Position',[70 70 1000 400]) %[0 0.547 0.841]
schedule=120;
plot(50:1:schedule+49,dataSet(1:schedule,1),'Color',[0.4 0.4 0.4],'LineWidth',2); hold on;
plot(50:1:schedule+49,dataSet(1:schedule,2),'Color',[0 0.547 0.841],'LineWidth',2); hold on;
plot(50:1:schedule+49,dataSet(1:schedule,3),'r','LineWidth',2); hold on;

schedule=50;

plot(1:1:schedule,dataSetTemp(1:schedule,1),'Color',[0.4 0.4 0.4],'LineWidth',2); hold on;
plot(1:1:schedule,dataSetTemp(1:schedule,2),'Color',[0 0.547 0.841],'LineWidth',2); hold on;
plot(1:1:schedule,dataSetTemp(1:schedule,3),'r','LineWidth',2); hold on;

ylabel('Métrica [%]');
xlabel('Iteración');
legend('Solución evaluada','Solución aceptada','Mejor solución hallada');
legend('Location','northeastoutside')
annotation('textbox',[.725 .66 .1 .1],'String',strcat('Máximo valor hallado:',{' '},num2str(max(dataSet(1:schedule,3)))));
set(gca,'FontSize',12)
set(gca,'GridAlpha',0.25)
set(gca,'Position',[0.10038 0.158 0.59843 0.797])
grid on
hold on
