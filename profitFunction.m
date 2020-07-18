
    join=zeros(0,2);
    ft2sm=zeros(0,2);

    for ii=1:nFT
        for jj=1:nSM
            
          facility = scenario.Children.Item(strcat('FT_',num2str(ii)));
          access = facility.GetAccess(strcat('Satellite/SM_',num2str(jj)));
          access.ComputeAccess();
          intervalCollection = access.ComputedAccessIntervalTimes;
          
          if intervalCollection.count ~= 0 % Existen contactos?
          join = sortrows(union(join,cell2mat(intervalCollection.ToArray(0, -1)),'rows'));
          end
        end
    end
    
    if size(join,1)>0
        tempVector = join(1,:);
            for k=2:size(join,1)
                if tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)>tempVector(2)
                   tempVector(2) = join(k,2);
                elseif tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)<=tempVector(2)
                   tempVector(2) = tempVector(2);
                else
                   ft2sm=union(ft2sm,tempVector,'rows');
                   tempVector = join(k,:);
                end
            end
        ft2sm=union(ft2sm,tempVector,'rows');
        M(1)=sum(ft2sm(:,2)-ft2sm(:,1));
        M(1)=(M(1)/scTime)*100;
    else
        M(1)=0;
    end

%% Profit w/relay

if nRelay >=1

%% Ground Assets -> Relay Net
    
    join=zeros(1,2);
    ft2relay=zeros(0,2);
    
    for ii=1:nFT
        for jj=1:nRelay
            
          facility = scenario.Children.Item(strcat('FT_',num2str(ii)));
          access = facility.GetAccess(strcat('Satellite/RELAY_',num2str(jj)));
          access.ComputeAccess();
          intervalCollection = access.ComputedAccessIntervalTimes;
          
              if intervalCollection.count ~= 0 % Existen contactos?
              join = sortrows(union(join,cell2mat(intervalCollection.ToArray(0, -1)),'rows'));
              end
        end
    end
    
        tempVector = join(1,:);
        for k=2:size(join,1)
            if  tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)>tempVector(2)
               tempVector(2) = join(k,2);
            elseif tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)<=tempVector(2)
               tempVector(2) = tempVector(2);  
            else
               ft2relay=union(ft2relay,tempVector,'rows');
               tempVector = join(k,:);
            end
        end
        ft2relay=union(ft2relay,tempVector,'rows');

%% Mission Satellites -> Relay net

    join=zeros(1,2);
    sm2relay=zeros(0,2);
    
    for ii=1:nSM
        for jj=1:nRelay
            
          sat = scenario.Children.Item(strcat('SM_',num2str(ii)));
          access = sat.GetAccess(strcat('Satellite/RELAY_',num2str(jj)));
          access.ComputeAccess();
          intervalCollection = access.ComputedAccessIntervalTimes;
          
              if intervalCollection.count ~= 0 % Existen contactos?
              join = sortrows(union(join,cell2mat(intervalCollection.ToArray(0, -1)),'rows'));
              end
        end
    end
    
        tempVector = join(1,:);
        for k=2:size(join,1)
            if  tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)>tempVector(2)
               tempVector(2) = join(k,2);
            elseif tempVector(1)<=join(k,1) && join(k,1)<=tempVector(2) && join(k,2)<=tempVector(2)
               tempVector(2) = tempVector(2);        
            else
               sm2relay=union(sm2relay,tempVector,'rows');
               tempVector = join(k,:);
            end
        end
        sm2relay=union(sm2relay,tempVector,'rows');

%% Intersections 

          % gs2relay INT ms2relay
          k=1; firstSet=ft2relay; secondSet=sm2relay; overlap=zeros(max([size(firstSet,1) size(secondSet,1)])-1,2);

            for ii=1:size(firstSet,1)
                for jj=1:size(secondSet,1)
                    if(firstSet(ii,1)<=secondSet(jj,1)&&firstSet(ii,2)>secondSet(jj,1)&&firstSet(ii,2)<=secondSet(jj,2))
                    overlap(k,:)=[secondSet(jj,1) firstSet(ii,2)]; k=k+1;
                    elseif(firstSet(ii,1)>=secondSet(jj,1)&&firstSet(ii,1)<secondSet(jj,2)&&firstSet(ii,2)>=secondSet(jj,2))
                    overlap(k,:)=[firstSet(ii,1) secondSet(jj,2)]; k=k+1;
                    elseif(firstSet(ii,1)>=secondSet(jj,1)&&firstSet(ii,2)<=secondSet(jj,2))
                    overlap(k,:)=[firstSet(ii,1) firstSet(ii,2)]; k=k+1;
                    elseif(firstSet(ii,1)<secondSet(jj,1)&&firstSet(ii,2)>secondSet(jj,2))
                    overlap(k,:)=[secondSet(jj,1) secondSet(jj,2)]; k=k+1;
                    end
                end
            end
          
          relay2msAndFt=overlap;
          
          % gs2relay INT gs2ms
          k=1; firstSet=ft2relay; secondSet=ft2sm; overlap=zeros(max([size(firstSet,1) size(secondSet,1)])-1,2);

            for ii=1:size(firstSet,1)
                for jj=1:size(secondSet,1)
                    if(firstSet(ii,1)<=secondSet(jj,1)&&firstSet(ii,2)>secondSet(jj,1)&&firstSet(ii,2)<=secondSet(jj,2))
                    overlap(k,:)=[secondSet(jj,1) firstSet(ii,2)]; k=k+1;
                    elseif(firstSet(ii,1)>=secondSet(jj,1)&&firstSet(ii,1)<secondSet(jj,2)&&firstSet(ii,2)>=secondSet(jj,2))
                    overlap(k,:)=[firstSet(ii,1) secondSet(jj,2)]; k=k+1;
                    elseif(firstSet(ii,1)>=secondSet(jj,1)&&firstSet(ii,2)<=secondSet(jj,2))
                    overlap(k,:)=[firstSet(ii,1) firstSet(ii,2)]; k=k+1;
                    elseif(firstSet(ii,1)<secondSet(jj,1)&&firstSet(ii,2)>secondSet(jj,2))
                    overlap(k,:)=[secondSet(jj,1) secondSet(jj,2)]; k=k+1;
                    end
                end
            end         

          ft2smAndRelay=overlap;
          
          M(1)=sum(sm2relay(:,2)-sm2relay(:,1))+sum(ft2relay(:,2)-ft2relay(:,1))+sum(ft2sm(:,2)-ft2sm(:,1))-sum(ft2smAndRelay(:,2)-ft2smAndRelay(:,1));
          M(1)=((M(1)/scTime)*100)/3;
          
end