%%  Metrica inicial

    profitFunction;
    
%%  Variables del algoritmo

    rng('shuffle');
    M = [0,0,0]; % [M_iteracion, M_aceptada, M_optima]
    dataSet=zeros(itMax,9,nSM+nRelay); % dataSet for plotting
    Ps(2,:,:)=Ps(1,:,:);
    Pft(2,:,:)=Pft(1,:,:);

%%  Configuraciones

    probabilities=zeros(itMax,2);
    delta = 0.01; % Delta esperado entre iteraciones
    pi = 0.75;   % Probability of accepting worse solution at the start
    pf = 0.001; % Probability of accepting worse solution at the end
    tmax = -delta/log(pi); % Initial temperature
    ti = tmax; % Temperatura de iteración
    tf = -delta/log(pf); % Final temperature
    alpha = (tf/ti)^(1/(itMax-1)); % Fractional reduction every cycle
    delta=0; itRst=0;
    element = 1; attempts = 0; nAsset = 1; itRst = 0;
    
%%  Algoritmo de búsqueda

for schedule=1:itMax

    if nAsset>nSM+nRelay+nFT
       nAsset=1; element=1;
    end
    if element>6
       element = 1; nAsset = nAsset+1;
    end
    while nAsset<=nSM+nRelay && dPs(nAsset,element)==0
        element = element + 1;
        if element>6
           element = 1; nAsset = nAsset+1;
        end
    end
    if nAsset > nSM+nRelay
        while nAsset<=nSM+nRelay+nFT && dFt(nAsset-nRelay-nSM)==0
            nAsset = nAsset+1;
            if nAsset>nSM+nRelay+nFT
            nAsset = 1; break;
            end
        end
        while nAsset<=nSM+nRelay && dPs(nAsset,element)==0
            element = element + 1;
            if element>6
               element = 1; nAsset = nAsset+1;
            end
        end
    end

    fprintf(strcat(num2str(schedule/itMax*100),' %')); fprintf('\n');
    fprintf(strcat('Optimizando parametro:_',num2str(element),'_del asset_', num2str(nAsset))); fprintf('\n');
    Ps(1,:,:) = Ps(2,:,:);
    
    %% Encontrar solución vecina

    if nAsset <= nSM + nRelay  % Finding neighbour in Flight Segment

    parametrosTemp = Ps(1,:,nAsset);
    parametrosTemp(element) = (A(element,2,nAsset)-A(element,1,nAsset))*rand + A(element,1,nAsset);

    keplerian.SizeShape.SemiMajorAxis = parametrosTemp(1);
    keplerian.SizeShape.Eccentricity = parametrosTemp(2);
    keplerian.Orientation.Inclination = parametrosTemp(3);
    keplerian.Orientation.AscNode.Value = parametrosTemp(4);
    keplerian.Orientation.ArgOfPerigee = parametrosTemp(5);
    keplerian.Location.Value = parametrosTemp(6);

    Ps(1,:,nAsset) = parametrosTemp;

    if nAsset<=nSM
        sat = scenario.Children.Item(strcat('SM_',num2str(nAsset)));
    else
        sat = scenario.Children.Item(strcat('RELAY_',num2str(nAsset-nSM))); 
    end

    sat.Propagator.InitialState.Representation.Assign(keplerian);
    sat.Propagator.Propagate;

    else   % Finding neighbour in Ground Segment
        
        parametrosTemp(1:4)=F(randi([1 size(F,1)],1),:);
        while ismember(parametrosTemp(1),Pft(2,1,:))% && dFt(nAsset-nSM-nRelay)>=0
        parametrosTemp(1:4)= F(randi([1 size(F,1)],1),:);
        end

        facility = scenario.Children.Item(strcat('FT_',num2str(nAsset-nSM-nRelay)));
        facility.Position.AssignGeodetic(parametrosTemp(2),parametrosTemp(3),parametrosTemp(4));
        Pft(1,:,nAsset-nSM-nRelay)=parametrosTemp(1:4);

    end
            
    %% Calculo de métrica
    
    profitFunction;
    delta=delta+abs(M(2)-M(1));
    
    %% Bloque de decisión
    
    if(M(1)>M(2)) % Solucion vecina mejor que solución actual?
        M(2)=M(1);
        Ps(2,:,:) = Ps(1,:,:);
        Pft(2,:,:) = Pft(1,:,:);
        if (M(1)>M(3)) % Solución vecina mejor que la solución optima?
            itRst = 0;
            M(3) = M(1);
            Ps(3,:,:) = Ps(1,:,:);
            Pft(3,:,:) = Pft(1,:,:);
        end
    else
        if ((exp((M(1)-M(2))/(ti)))>rand()) % Criterio de acept.
            probabilities(schedule,:)=[1,(exp((M(1)-M(2))/(ti)))];
            M(2)=M(1);
            Ps(2,:,:) = Ps(1,:,:);
            Pft(2,:,:) = Pft(1,:,:);
        else % Rechazo el movimiento
            probabilities(schedule,:)=[2,(exp((M(1)-M(2))/(ti)))];
        end
    end

    ti=alpha*ti;
    
    if(M(1)==M(2))
        M(2)=M(1);
        Ps(2,:,:) = Ps(1,:,:);
        Pft(2,:,:) = Pft(1,:,:);
    else
        if (nAsset <= (nSM+nRelay))     % Finding neighbour in Flight Segment
            element = element + 1;
        else                            % Finding neighbour in Ground Segment
            nAsset = nAsset+1;
        end
        if nAsset>nSM+nRelay+nFT
            nAsset=1;
        end
    end
    
    for k=1:nSM+nRelay
        dataSet(schedule,1:3,k)=[M(1) M(2) M(3)]; dataSet(schedule,4:end,k)=Ps(1,:,k);
    end
    
    itRst = itRst+1;
    if itRst>rstAt
        Ps(2,:,:) = Ps(3,:,:);
        Pft(2,:,:) = Pft(3,:,:);       
        itRst=0;
        fprintf('RESET');
    end
    
    %averageTime(schedule,1)=toc;
    
end