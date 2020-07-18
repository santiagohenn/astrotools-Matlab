%%  Metrica inicial

    profitFunction;
    
%%  Variables del algoritmo

    rng('shuffle');
    M = [0,0,0]; % [M_iteracion, M_aceptada, M_optima]
    dataSet=zeros(itMax,9,nSM+nRelay); % dataSet for plotting
    Ps(2,:,:)=Ps(1,:,:); Ps(3,:,:)=Ps(1,:,:);
    Pft(2,:,:)=Pft(1,:,:); Pft(3,:,:)=Pft(1,:,:);
    gradient = 1; % 0->left 1->right
    attempts = 0; stuck = 0; decStep = 0; maxStuck = 3;
    
%%  Which asset and element?
    
    nAsset = 1; element = 6; 

%%  Schedule
    
    itMax = round((A(element,2,nAsset)-A(element,1,nAsset)) / dPs(nAsset,element))
    
for schedule=1:itMax
    
    fprintf(strcat(num2str(schedule/itMax*100),' %')); fprintf('\n');
    Ps(1,:,:) = Ps(2,:,:);
    
    %% Encontrar solución vecina
    
    parametrosTemp = Ps(1,:,nAsset);
    parametrosTemp(element)=parametrosTemp(element)+gradient*dPs(nAsset,element)

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
            
    %% Profit Calculation
    profitFunction;
    
    if(M(1)>M(2))    % Neighbour Solution Better?
        M(2)=M(1); M(3)=M(1);
        Ps(3,:,:) = Ps(2,:,:);
        Pft(3,:,:) = Pft(2,:,:);
    end
    
    % No matter what, i keep evaluating:
    Ps(2,:,:) = Ps(1,:,:);
    Pft(2,:,:) = Pft(1,:,:); 
    
    for k=1:nSM+nRelay
        dataSet(schedule,1:3,k)=[M(1) M(2) M(3)]; dataSet(schedule,4:end,k)=Ps(1,:,k);
    end

end

