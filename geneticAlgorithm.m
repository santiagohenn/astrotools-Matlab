%% Variables del algoritmo

rng('shuffle');
population=80;
discarted=population/2;
genes=6;
generation=1; 
generations=itMax;
pC=0.95;               % probabilidad de crossover
pM=0.2;               % probabilidad de mutación
crossOverPoint=3;
phenotypes=zeros(population,genes,nSM+nRelay+nFT);  %[P1, P2, P3, ... , P_Ngenes]
allele=zeros(1,6);
fitness=zeros(population,1);
totalFitness=zeros(generations,1);
maxFitness=zeros(generations,1);
bestIndividual=zeros(1,genes+1,nSM+nRelay+nFT);

    %% Initial population

        for k=1:nSM+nRelay
            for i=1:population
                for j=1:6
                    phenotypes(i,j,k)=(A(j,2,k)-A(j,1,k))*rand+A(j,1,k);
                end
            end
        end

        for k=nSM+nRelay+1:nSM+nRelay+nFT
            for i=1:population
                if ~(dFt(k-nSM-nRelay) == 0)
                    phenotypes(i,1:4,k)=F(randi([1 size(F,1)],1),:);
                else
                    phenotypes(i,1:4,k)=F(k-nSM-nRelay,:);
                end
            end
        end

%% Fitness Computation

for generation=1:generations
    
    fprintf(strcat(num2str(generation/generations*100),' %')); fprintf('\n');
    
    if generation==1
        from = 1;
    else
        from = population-discarted;
    end
    
for individuo=from:population

    fprintf(strcat(num2str(individuo/population*100),' %')); 

    for k=1:nSM+nRelay
        if ~(1 && all(dPs(k,:) == 0)) % Check if algun parametro es variable
            allele=phenotypes(individuo,:,k);
            keplerian.SizeShape.SemiMajorAxis = allele(1);
            keplerian.SizeShape.Eccentricity = allele(2);
            keplerian.Orientation.Inclination = allele(3);
            keplerian.Orientation.AscNode.Value = allele(4);
            keplerian.Orientation.ArgOfPerigee = allele(5);
            keplerian.Location.Value = allele(6);

            if k<=nSM
                sat = scenario.Children.Item(strcat('SM_',num2str(k)));
            else
                sat = scenario.Children.Item(strcat('RELAY_',num2str(k-nSM))); 
            end        

            sat.Propagator.InitialState.Representation.Assign(keplerian);
            sat.Propagator.Propagate;
        end
    end

    for k=1:nFT
        allele=phenotypes(individuo,:,nSM+nRelay+k);
        facility = scenario.Children.Item(strcat('FT_',num2str(k)));
        facility.Position.AssignGeodetic(allele(2),allele(3),allele(4));
    end
    
    profitFunction;
    fitness(individuo,1)=M(1);
%     fprintf(strcat('\n','metrica_',num2str(fitness(individuo,1)),'\n'));
end

%% Selection (by fitness proportionate selection)

    totalFitness(generation)=sum(fitness(:,1));
    spec=zeros(population,genes+1,nSM+nRelay+nFT);
    avgFitness=totalFitness/population;
    maxFitness(generation)=max(fitness);

    for k=1:nSM+nRelay+nFT
        spec(:,:,k)=sortrows(horzcat(fitness/totalFitness(generation),phenotypes(:,:,k)),-1);
    end
    
    fitness(:,1)=spec(:,1,1)*totalFitness(generation);
    
    if (spec(1,1,1)>bestIndividual(1,1,1))
        bestIndividual(1,:,:)=spec(1,:,:);
    end
    
    if generation==generations
        break;
    end

%% Normalizing

    for k=1:nSM+nRelay+nFT
        for j=2:population
            spec(j,1,k)=spec(j-1,1,k)+spec(j,1,k);
        end
    end

%% Selecting

selection=zeros(population,genes+1,nSM+nRelay+nFT);
for j=1:population-discarted
    chance=rand;
    for indiv=1:population
        if (spec(indiv,1,1)>=chance)          % CONTRAOLAR ESTOOOOOO
        selection(j,:,:)=spec(indiv,:,:);     % Using UNION avoids duplicates
        break
        end
    end
end

%%
%     for k=1:nSM
%         for i=1:population
%             accumFitness(i,1,k)=sum(spec(1:i,1,k)); % Accumulated normalized fitness values
%         end
%     end

%%
%     for k=1:nSM
%         for i=1:population-discarted
%             prob=rand;
%             for j=1:population
%                 if (accumFitness(j,1)>prob)
%                     selection(i,:,k)=accumFitness(j,:,k);     % Using UNION avoids duplicates
%                     break;
%                 end
%             end    
%         end
%     end

%% Offspring - Crossover

for k=1:nSM+nRelay

    for individuo=1:(population-discarted)/2
        offspring=zeros(2,genes+1);

        if(rand<pC)
            crossOverPoint=randi([2 genes],1);
            offspring(1,1:crossOverPoint)=selection(2*individuo,1:crossOverPoint,k);
            offspring(2,1:crossOverPoint)=selection(2*individuo-1,1:crossOverPoint,k);
            offspring(1,crossOverPoint+1:end)=selection(2*individuo-1,crossOverPoint+1:end,k);
            offspring(2,crossOverPoint+1:end)=selection(2*individuo-1,crossOverPoint+1:end,k);
            
%             for j=2:crossOverPoint
%                 offspring(1,j)=selection(2*individuo,j,k);
%                 offspring(2,j)=selection(2*individuo-1,j,k);
%             end
%             for j=crossOverPoint:genes+1
%                 offspring(1,j)=selection(2*individuo-1,j,k);
%                 offspring(2,j)=selection(2*individuo,j,k);
%             end

        else
            offspring(1,:)=selection(2*individuo-1,:,k);
            offspring(2,:)=selection(2*individuo,:,k);            
        end

        selection((population-discarted+2*individuo),:,k)=offspring(1,:);
        selection((population-discarted+2*individuo-1),:,k)=offspring(2,:);
        
    end
end

    for individuo=1:(population-discarted)/2
        offspringFT=zeros(2,7,nSM+nRelay+nFT);
        if(rand<pC)
            crossOverPoint=randi([nSM+nRelay+1 nSM+nRelay+nFT],1);
            offspringFT(1,:,1:crossOverPoint)=selection(2*individuo,:,1:crossOverPoint);
            offspringFT(1,:,crossOverPoint:end)=selection(2*individuo-1,:,crossOverPoint:end);
            offspringFT(2,:,1:crossOverPoint)=selection(2*individuo,:,1:crossOverPoint);
            offspringFT(2,:,crossOverPoint:end)=selection(2*individuo-1,:,crossOverPoint:end);
        else
            offspringFT(1,:,:)=selection(2*individuo-1,:,:);
            offspringFT(2,:,:)=selection(2*individuo,:,:);            
        end
            selection((population-discarted+2*i),:,:)=offspringFT(1,:,:);
            selection((population-discarted+2*i-1),:,:)=offspringFT(2,:,:);
    end

%% Offspring - Mutation

    for k=1:nSM+nRelay
        for i=population-discarted:population
            for j=1:genes
                if(rand<pM)
                    selection(i,j+1,k)=(A(j,2,k)-A(j,1,k))*rand+A(j,1,k);
                end
            end
        end   
        phenotypes(i,:,k)=selection(i,2:end,k);
    end

    for k=1:nFT
        if (dFt(k) > 0)
            for i=population-discarted:population
                if(rand<pM)
                    allele(1:4)=F(randi([1 size(F,1)],1),:);
                    while ismember(allele(1),selection(i,1,nSM+nRelay:end))% && dFt(nAsset-nSM-nRelay)>=0
                    allele(1:4)= F(randi([1 size(F,1)],1),:);
                    end
                phenotypes(i,:,nSM+nRelay+k)=allele;
                end
            end
        end
    end

end

%% Final Fitness Computation
%     
% for i=1:population
%     
%     fprintf(strcat(num2str(i/population*100),' %')); fprintf('\n');
%     for k=1:nSM
%         allele=phenotypes(i,:,k);
%         sat = scenario.Children.Item(strcat('SM_',num2str(k)));   
%         keplerian.Orientation.Inclination = allele(3);
%         keplerian.Orientation.AscNode.Value = allele(4);
%         keplerian.Orientation.ArgOfPerigee = allele(6);
%         keplerian.Location.Value = 0;
%         sat.Propagator.InitialState.Representation.Assign(keplerian);
%         sat.Propagator.Propagate;
%     end
%     
%     profitSimple;
%     M(1);
%     
% end
% 
% %     totalFitness(generation)=sum(fitness(:,1));
% %     spec=zeros(population,genes+1,nSM);
% %     avgFitness=totalFitness/population;
% %     maxFitness=max(fitness);
% 
%     for k=1:nSM
%         spec(:,:,k)=sortrows(horzcat(fitness/totalFitness(generation),phenotypes(:,:,k)),-1);
%     end  
%     
%%

    Ps(3,:,1:nSM+nRelay) = bestIndividual(1,2:end,1:nSM+nRelay);
    Pft(3,1:4,1:nFT) = bestIndividual(1,2:5,nSM+nRelay+1:end);
    profitFunction;
    M(3)=max(maxFitness);
    fprintf('\n'); fprintf('DONE'); fprintf('\n');
    
    %%

