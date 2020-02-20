clear

folder_content = dir("Connectomes/");
control_index = 1;
pd_icd_index = 1;
pd_no_icd_index = 1;

%setti a 0 tutte le connessioni che non sopravvivono alla treshold
for i = 1:length(folder_content)
    if startsWith(folder_content(i).name, "HC_")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        control(:,:,control_index) = temp + temp' - diag(diag(temp));
        control_index = control_index + 1;
    elseif startsWith(folder_content(i).name, "PD_ICD")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        pd_icd(:,:,pd_icd_index) = temp + temp' - diag(diag(temp));
        pd_icd_index = pd_icd_index + 1;
    elseif startsWith(folder_content(i).name, "PD_NO_ICD")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        pd_no_icd(:,:,pd_no_icd_index) = temp + temp' - diag(diag(temp));
        pd_no_icd_index = pd_no_icd_index + 1;
    end
end

%figure, imagesc(control(:,:,1)), axis equal tight, colorbar;
%figure, imagesc(pd_icd(:,:,1)), axis equal tight, colorbar;
%figure, imagesc(pd_no_icd(:,:,1)), axis equal tight, colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge level statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%

th = 0;
nc = size(control,3);
group_presence = sum(control > th,3) / nc; %ritorna la percentuale di soggetti che hanno la connessione

for i = 1:size(group_presence,1)
    for j = 1:size(group_presence,2)
        if group_presence(i,j) < 0.7
            control(i,j,:) = 0;
        end
    end
end


npi = size(pd_icd,3);
group_presence = sum(pd_icd > th,3) / npi;

for i = 1:size(group_presence,1)
    for j = 1:size(group_presence,2)
        if group_presence(i,j) < 0.7
            pd_icd(i,j,:) = 0;
        end
    end
end


npni = size(pd_no_icd,3);
group_presence = sum(pd_no_icd > th,3) / npni;

for i = 1:size(group_presence,1)
    for j = 1:size(group_presence,2)
        if group_presence(i,j) < 0.7
            pd_no_icd(i,j,:) = 0;
        end
    end
end

control_pd_icd_stat_differences = [];
control_pd_no_icd_stat_differences = [];
% Differences between (controls and pd_icd) and (controls and pd_no_icd)
for i = 1:size(control,1)
    for j = 1:size(control,2)
        % Test only the lower triangular
        if j <= i
            a = squeeze(control(i,j,:)); %array dei valori di tutti i pazienti della connessione dal nodo i al nodo j
            b = squeeze(pd_icd(i,j,:));
            c = squeeze(pd_no_icd(i,j,:));
            % is the connection in pc_icd weaker than in control?
            [H1,P1] = ttest2(a,b,'tail','right');%controllo più forte dei casi? P1 pvalue
            if H1 == 1 % se rifiuto H0 (che è i casi sono più forti o uguali dei controlli)
                control_pd_icd_stat_differences = [control_pd_icd_stat_differences;[i,j,P1]];
            end
            % is the connection in pc_no_icd weaker than in control?
            [H2,P2] = ttest2(a,c,'tail','right');
            if H2 == 1
                control_pd_no_icd_stat_differences = [control_pd_no_icd_stat_differences;[i,j,P2]];
            end
        end
    end
end

% Number of comparisons
bonf_coefficient = ((size(control,1)*size(control,1))-size(control,1))/2 + size(control,1); %quantità di test che ho fatto, mi serve per bonferroni
control_pd_icd_stat_differences_bonf = [];
control_pd_no_icd_stat_differences_bonf = [];
for i = 1:size(control_pd_icd_stat_differences,1)
    if control_pd_icd_stat_differences(i,3) < 0.05/bonf_coefficient
        control_pd_icd_stat_differences_bonf = [control_pd_icd_stat_differences_bonf;[control_pd_icd_stat_differences(i,:)]];
    end
end

for i = 1:size(control_pd_no_icd_stat_differences,1)
    if control_pd_no_icd_stat_differences(i,3) < 0.05/bonf_coefficient
        control_pd_no_icd_stat_differences_bonf = [control_pd_no_icd_stat_differences_bonf;[control_pd_no_icd_stat_differences(i,:)]];
    end
end

% None of the edges seems to survive bonferroni correction :(



%%%%%%%%%%%
% NBS
%%%%%%%%%%%
UI.method.ui='Run NBS'; 
UI.test.ui='t-test';
UI.size.ui='Extent';
UI.thresh.ui='3.1';
UI.perms.ui='5000';
UI.alpha.ui='0.05';
UI.contrast.ui='[-1,1]'; 
%UI.design.ui='SchizophreniaExample/designMatrix.txt';
UI.exchange.ui=''; 
UI.matrices.ui=control_pd_icd;
%UI.node_coor.ui='SchizophreniaExample/COG.txt';                         
%UI.node_label.ui='SchizophreniaExample/nodeLabels.txt';

%NBSrun(UI,[])...
%global nbs
%explore nbs structure..
