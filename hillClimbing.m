%%  Metrica inicial

    profitFunction;
    
%%  Variables del algoritmo

    rng('shuffle');
    M = [0,0,0]; % [M_iteracion, M_aceptada, M_optima]
    dataSet=zeros(itMax,9,nSM+nRelay); % dataSet for plotting
    Ps(2,:,:)=Ps(1,:,:); Ps(3,:,:)=Ps(1,:,:);
    Pft(2,:,:)=Pft(1,:,:); Pft(3,:,:)=Pft(1,:,:);
    gradient = 1; % 0->left 1->right
    element = 1; attempts = 0; nAsset = 1; stuck = 0; decStep = 0; maxStuck = 3;

%%  Algoritmo de búsqueda
    
    while dPs(nAsset,element)==0 && nAsset<=nSM+nRelay
        element = element + 1;
        if element>6
           element = 1; nAsset = nAsset+1;
        end
    end
    if nAsset > nSM+nRelay
        while dFt(nAsset-nRelay-nSM)==0
            nAsset = nAsset+1;
            if nAsset>nSM+nRelay+nFT
               itMax=0; break;
            end
        end
    end
    
for schedule=1:itMax
    
    fprintf(strcat(num2str(schedule/itMax*100),' %')); fprintf('\n');
    fprintf(strcat('Optimizando parametro:_',num2str(element),'_del asset_', num2str(nAsset))); fprintf('\n');
    Ps(1,:,:) = Ps(2,:,:);
    
    %% Encontrar solución vecina

    if nAsset <= nSM + nRelay  % Finding neighbour in Flight Segment

    parametrosTemp = Ps(1,:,nAsset);
    parametrosTemp(element)=parametrosTemp(element)+gradient*dPs(nAsset,element);

        if element<=3 || (A(element,2,nAsset)-A(element,1,nAsset))~=360
            if parametrosTemp(1,element)<A(element,1,nAsset)
               parametrosTemp(1,element)=A(element,1,nAsset);
            elseif parametrosTemp(1,element)>A(element,2,nAsset)
               parametrosTemp(1,element)=A(element,2,nAsset);
            end
        else
            if parametrosTemp(1,element) < 0
               parametrosTemp(1,element) = 360 + parametrosTemp(1,element);
            elseif parametrosTemp(1,element) > 360
               parametrosTemp(1,element) = parametrosTemp(1,element) - 360;
            end
        end
        
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
    parametrosTemp(1:4) = Pft(2,:,nAsset-nSM-nRelay);

        while ismember(parametrosTemp(1),Pft(1,1,:)) && dFt(nAsset-nSM-nRelay)~=0
            if (parametrosTemp(1)+gradient)>size(F,1)
                parametrosTemp(1:4) = F(1,:);
            elseif (parametrosTemp(1)+gradient)<1
                parametrosTemp(1:4) = F(size(F,1),:);
            else
                parametrosTemp(1:4) = F(parametrosTemp(1)+gradient,:);
            end
        end

    facility = scenario.Children.Item(strcat('FT_',num2str(nAsset-nSM-nRelay)));
    facility.Position.AssignGeodetic(parametrosTemp(2),parametrosTemp(3),parametrosTemp(4));
    Pft(2,:,nAsset-nSM-nRelay)=parametrosTemp(1:4);

    end
            
    %% Profit Calculation
    profitFunction;
    
    %% Decision Block (Hill Climbing)
    
    if(M(1)>M(2))    % Neighbour Solution Better?
        M(2)=M(1); M(3)=M(1);
        Ps(2,:,:) = Ps(1,:,:); Ps(3,:,:) = Ps(2,:,:);
        Pft(2,:,:) = Pft(1,:,:); Pft(3,:,:) = Pft(2,:,:);
    else
        
    if attempts < 1
        gradient=(-1)*gradient;
        attempts = attempts + 1;
    else
        attempts = 0; element = element + 1; decStep = decStep + 1;
        if element>6
            nAsset=nAsset+1; element=1;
        elseif nAsset<=nSM+nRelay
            while (nAsset<=nSM+nRelay && dPs(nAsset,element)==0)
                element = element+1;
                if element>6
                    element = 1; nAsset = nAsset+1;
                end
            end
        elseif nAsset>nSM+nRelay
            while (nAsset>nSM+nRelay && dFt(nAsset-nRelay-nSM)==0)
                nAsset = nAsset+1;
                if nAsset>nSM+nRelay+nFT
                   nAsset = 1;
                end
            end
            while (nAsset<=nSM+nRelay && dPs(nAsset,element)==0)
                element = element+1;
                if element>6
                    element = 1; nAsset = nAsset+1;
                end
            end
        end
    end
    
    end
    
    for k=1:nSM+nRelay
        dataSet(schedule,1:3,k)=[M(1) M(2) M(3)]; dataSet(schedule,4:end,k)=Ps(1,:,k);
    end
    
    if decStep>=(nSM+nRelay)*6+nFT
        fprintf('\n'); fprintf('STEP DECREASED'); fprintf('\n');
        dPs = dPs/2; decStep = 0;
        stuck = stuck + 1;
    else 

    end
    if stuck==3
        fprintf('\n'); fprintf('Local Optimum Reached_');
        break;
    end

end

