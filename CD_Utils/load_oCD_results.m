function [t,n,ib,Cplot,Cplotall,Zvg,Qvg,Ncvg] = ...
                load_oCD_results(n,subjects,missing_subjs,t,p,ts,...
                                     goRange,tag,inpath)
    if missing_subjs == [];
        ib = loadib(subjects);
    else
        ib = loadib(subjects,missing_subjs);
    end
       
    krange = 1:numel(ib); % index over people
    Zvg = zeros(size(goRange,1),numel(ib),4,t);  % Z, Zall
    Qvg = zeros(size(goRange,1),numel(ib),5,t);  % Q, Qall, Qconsall
    Ncvg = zeros(size(goRange,1),numel(ib),4,t+1); % Nc, Ncall
    Cplot = zeros(size(goRange,1),p*numel(ib),n,t); %
    Cplotall = zeros(size(goRange,1),p*numel(ib),n,t);

for goix=1:size(goRange,1);
    gamma = goRange(goix,1);
    omega = goRange(goix,2);

for k=krange
    person = ib(k);
    Trange = 1:t;
    disp(['SUBJECT ' num2str(person)]);
    
    % read in community detection results    
    Call = csvread([inpath '/subj' num2str(person) '_' num2str(gamma) '_' ...
            num2str(omega) '_ts' num2str(ts) tag '_Call.txt']);
    Cnewall = csvread([inpath '/subj' num2str(person) '_' num2str(gamma) '_' ...
            num2str(omega) '_ts' num2str(ts) tag '_Cnewall.txt']);
    Qs = csvread([inpath '/subj' num2str(person) '_' num2str(gamma) '_' ...
            num2str(omega) '_ts' num2str(ts) tag '_Qs.txt']);
    Zs = csvread([inpath '/subj' num2str(person) '_' num2str(gamma) '_' ...
            num2str(omega) '_ts' num2str(ts) tag '_Zs.txt']);
    Qconsall = csvread([inpath '/subj' num2str(person) '_' num2str(gamma) '_' ...
            num2str(omega) '_ts' num2str(ts) tag '_Qconsall.txt']);
        
        %[size1,size2] = size(Call);
        %assert(size2==n);
        %p = size1/t;
        %assert(mod(p,1)==0);
            for T = Trange
                if omega~=0 % for multislice reading of Zs and Qs
                % save info for each person and gammapair and slice (runwise slices, non-runwise Zs and Qs)
                Zvg(goix,k,:,1) = shiftdim([mean(Zs,1),var(Zs,0,1)],-1);  % Z, Zall
                Qvg(goix,k,:,1) = shiftdim([mean(Qs,1),var(Qs,0,1),Qconsall(1)],-1);  % Q, Qall, Qconsall
                else
                % save info for each person and gammapair and slice (one or runwise slices)
                Zvg(goix,k,:,T) = shiftdim([mean(Zs((T-1)*p*(p-1)/2+1:T*p*(p-1)/2,:),1),...
                                    var(Zs((T-1)*p*(p-1)/2+1:T*p*(p-1)/2,:),0,1)],-1);  % Z, Zall
                Qvg(goix,k,:,T) = shiftdim([mean(Qs((T-1)*p+1:T*p,:),1),var(Qs((T-1)*p+1:T*p,:),0,1),...
                                    Qconsall(T)],-1);  % Q, Qall, Qconsall
                end
                Ncvg(goix,k,:,T) = shiftdim([mean(max(Call((T-1)*p+1:T*p,:),[],2),1),...
                                    mean(max(Cnewall((T-1)*p+1:T*p,:),[],2),1),...
                                    var(max(Call((T-1)*p+1:T*p,:),[],2),0,1),...
                                    var(max(Cnewall((T-1)*p+1:T*p,:),[],2),0,1)],-1); % Nc, Ncall

                Cplot(goix,p*(k-1)+1:p*k,:,T) = shiftdim(Call((T-1)*p+1:T*p,:),-1);
                Cplotall(goix,p*(k-1)+1:p*k,:,T) = shiftdim(Cnewall((T-1)*p+1:T*p,:),-1);
            end % end loop over slices
            C3D = zeros(p,n,t);
            Cnew3D = zeros(p,n,t);
            for T=1:t
                C3D(:,:,T) = Call((T-1)*p+1:T*p,:);
                Cnew3D(:,:,T) = Cnewall((T-1)*p+1:T*p,:);
            end
            Ncvg(goix,k,:,Trange(end)+1) = shiftdim([mean(max(max(C3D,[],3),[],2),1),...
                                mean(max(max(Cnew3D,[],3),[],2),1),...
                                var(max(max(C3D,[],3),[],2),0,1),...
                                var(max(max(Cnew3D,[],3),[],2),0,1)],-1); % Nc overall
end % end loop over people

end % end loop over gammapairs

end