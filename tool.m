%% Variables generales
    format shortG;
    RT = 6378.137;       % Radio Terrestre (máximo) [Km]
    uT = 3.986004418E14; % Earth's Standard Gravitational parameter

%% Configuraciones generales
    
    itMax=100; rstAt=500;

%% Configuración del escenario

    root.UnitPreferences.Item('DateFormat').SetCurrentUnit('UTCG');
    root.UnitPreferences.SetCurrentUnit('Distance','Km');
    %scenario.SetTimePeriod('20 Jun 2020 16:00:05.000','20 Dec 2020 16:00:05.000'); root.Rewind;
    scenario.SetTimePeriod('20 Jun 2021 16:00:00.000','20 Dec 2021 16:00:00.000'); %root.Rewind;
    root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); % .SetCurrentUnit('JDate');
    root.Rewind;
    scTime=scenario.StopTime;
    root.ExecuteCommand('SubObjUnload *');
    
%% Espacios solución

    % Satélites
    
     % A(:,:,1) = [RT+400 RT+600; 0 0; 0 180; 0 360; 0 0; 0 360];
    A(:,:,1) = [7007.7 7007.7; 0 0; 97.8 97.8; 135 160; 0 0; 0 360];

    % Sun Synchr. mision:
    % A(:,:,1) = [7007.7 7007.7; 0 0; 97.9 97.9; 193 193; 90 90; 0 0];

    % Delta de paso para satélites
    dPs(1,:) = [0 0 0 0.2 0 0.2];
%     dPs(2,:) = [0 0 0 10 0 0];
%     dPs(2,:) = [0 0 0 0 0 0];
%     dPs(3,:) = [0 0 0 0 0 0];
%     dPs(4,:) = [0 0 5 5 0 0];
    
    % Facilidades en Tierra
    F=csvread('panamaCanal.csv');
    
    % Facilidades variables
    dFt = zeros(1,9);
    % dFt = [0 0 0 0 0 0 0 0 0];
    % dFt = [1 1 1 1 1 1];
    
%% Soluciones iniciales
    
    % Satélites
    nSM=1; nRelay=0; nFT=9; Ps=zeros(3,6,nSM+nRelay); Pft=zeros(3,4,nFT);
    
    % Misión
    % Ps(1,:,1)=[7007.7; 0; 97.8; 150; 0; 240];
    
    for i=1:nSM+nRelay
        for element=1:6
        Ps(1,element,i) = (A(element,2,i)-A(element,1,i))*rand + A(element,1,i);
        end
    end
    
%         for j=1:nFT
%         parametrosTemp(1:4)=F(randi([1 size(F,1)],1),:);
%         while ismember(parametrosTemp(1:4),Pft(1,1:4,:))
%         parametrosTemp(1:4)= F(randi([1 size(F,1)],1),:);
%         end
%         Pft(1,1:4,j)=parametrosTemp(1:4);
%     end
    
%    Facilidades
    for j=1:nFT
        Pft(1,:,j)=F(j,:);
    end
    
%     Pft(1,:,1)=F(randi([1 size(F,1)],1),:);
%     for j=2:nFT
%         parametrosTemp(1:4) = F(randi([1 size(F,1)],1),:);
%         while ismember(parametrosTemp(1),Pft(1,1,:)) && dFt(j)>=0
%         parametrosTemp(1:4)= F(randi([1 size(F,1)],1),:);
%         end
%         Pft(1,:,j)=parametrosTemp(1:4);
%     end
%     
    % Relay
    % Ps(1,:,3)=[7378.137,0,45,193.4537,180,0];
    % Ps(1,:,4)=[0,1000,80,0,0,0];
    % Ps(1,:,5)=[0,1000,80,0,0,0];

    % Facilidades en Tierra
    % Pft(1,:,1)=[1 -31.5263 -64.4674 0];
    % Pft(1,:,2)=[2 -54.5104 -67.1164 0];
    % Pft(1,:,3)=[3 49.9099 12.3955 0];
    
    % Pft(1,:,1)=[1 0 0 1];
    % Pft(1,:,2)=[2 90 0 1];
    % Pft(1,:,2)=[3 -90 0 1];
    % Pft(1,:,2)=[4 0 -180 1];
    % Pft(1,:,2)=[5 0 180 1];
    % Pft(1,:,2)=[6 0 180 1];
    
    % Pft(1,:,3)=[3 49.9099 12.3955 0];
    % Pft(1,:,4)=[4 49.9099 12.3955 0];
        
%     Pft(1,:,4)=F(randi([1 size(F,1)],1),:);
%     Pft(1,:,5)=F(randi([1 size(F,1)],1),:);
%     Pft(1,:,6)=F(randi([1 size(F,1)],1),:);

%% Crear elementos de la misión en el escenario
    
    fprintf('Loading Assets..'); fprintf('\n');

for nAsset=1:nSM+nRelay
    parametrosTemp=Ps(1,:,nAsset);
    if nAsset<=nSM
        sat = scenario.Children.New('eSatellite',strcat('SM_',num2str(nAsset)));
    else
        sat = scenario.Children.New('eSatellite',strcat('RELAY_',num2str(nAsset-nSM)));
    end
    keplerian = sat.Propagator.InitialState.Representation.ConvertTo('eOrbitStateClassical');
    sat.Propagator.step=5; % 60 segundos step 
    %keplerian.SizeShapeType = 'eSizeShapeAltitude';
    keplerian.SizeShapeType = 'eSizeShapeSemimajorAxis';
    keplerian.LocationType = 'eLocationTrueAnomaly';
    %keplerian.LocationType = 'eLocationMeanAnomaly';
    if ((parametrosTemp(1)<=42164) && (parametrosTemp(1)>=42162))
        keplerian.Orientation.AscNodeType = 0;  % 0=LAN
        fprintf('_location_type_LAN_');
    else
        keplerian.Orientation.AscNodeType = 1;  % 1=RAAN
        fprintf('_location_type_RAAN_');
    end
    keplerian.SizeShape.SemiMajorAxis = parametrosTemp(1);
    keplerian.SizeShape.Eccentricity = parametrosTemp(2);
    keplerian.Orientation.Inclination = parametrosTemp(3);
    keplerian.Orientation.AscNode.Value = parametrosTemp(4);
    keplerian.Orientation.ArgOfPerigee = parametrosTemp(5);
    keplerian.Location.Value = parametrosTemp(6);
    sat.Propagator.InitialState.Representation.Assign(keplerian);
    sat.SetPropagatorType(1);   % 1 = J2 perturbation, 2 = J4, 7 = Two body.
    sat.Propagator.Propagate;
    accessConstraints = sat.AccessConstraints;
    minAngle = accessConstraints.AddConstraint('eCstrGroundElevAngle');
    minAngle.EnableMin = true;
    minAngle.Min = 5; % Grados
    minGrazingAlt = accessConstraints.AddConstraint('eCstrGrazingAlt');
    minGrazingAlt.EnableMin = true;
    minGrazingAlt.Min = 100; % Km
end

for j=1:nFT
    facility = scenario.Children.New('eFacility',strcat('FT_',num2str(j)));
    facility.Position.AssignGeodetic(Pft(1,2,j),Pft(1,3,j),Pft(1,4,j));
end
    
    fprintf('..Assets Loaded.'); fprintf('\n');

%% Algoritmo de búsqueda
    
    simulatedAnnealing
    %hillClimbing
    %busquedaAleatoria
    %geneticAlgorithm
    %bruteForce
    
%% Presento la mejor solución

for nAsset=1:nSM+nRelay
    parametrosTemp = Ps(3,:,nAsset);
    if nAsset<=nSM
        sat = scenario.Children.Item(strcat('SM_',num2str(nAsset)));
    else
        sat = scenario.Children.Item(strcat('RELAY_',num2str(nAsset-nSM)));
    end
    keplerian.SizeShape.SemiMajorAxis = parametrosTemp(1);
    keplerian.SizeShape.Eccentricity = parametrosTemp(2);
    keplerian.Orientation.Inclination = parametrosTemp(3);
    keplerian.Orientation.AscNode.Value = parametrosTemp(4);
    keplerian.Orientation.ArgOfPerigee = parametrosTemp(5);
    keplerian.Location.Value = parametrosTemp(6);
    sat.Propagator.InitialState.Representation.Assign(keplerian);
    sat.Propagator.Propagate;
end

for nAsset=1:nFT
    parametrosTemp(1:4)=Pft(3,:,nAsset);
    facility = scenario.Children.Item(strcat('FT_',num2str(nAsset)));
    facility.Position.AssignGeodetic(parametrosTemp(2),parametrosTemp(3),parametrosTemp(4));
end   
    profitFunction;
    fprintf(strcat('Mejor valor encontrado_', num2str(M(3))));
    
%%  Refino
    
    %Ps(1,:,:) = Ps(3,:,:);
    %Ps(2,:,:) = Ps(3,:,:);
    %hillClimbing;

%%
