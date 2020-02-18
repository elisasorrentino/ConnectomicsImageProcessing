%cd into connectomes root directory

folder_content = dir;
control_index = 1;
pd_icd_index = 1;
pd_no_icd_index = 1;
th = 1;

for i = 1:length(folder_content)
    if startsWith(folder_content(i).name, "HC_")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        control(:,:,control_index) = temp + temp' - diag(diag(temp)); %#ok<SAGROW>
        control_index = control_index + 1;
    elseif startsWith(folder_content(i).name, "PD_ICD")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        pd_icd(:,:,pd_icd_index) = temp + temp' - diag(diag(temp)); %#ok<SAGROW>
        pd_icd_index = pd_icd_index + 1;
    elseif startsWith(folder_content(i).name, "PD_NO_ICD")
        temp = load(folder_content(i).folder + "/" + folder_content(i).name + "/Connectome/" + folder_content(i).name + "_connectome.csv");
        pd_no_icd(:,:,pd_no_icd_index) = temp + temp' - diag(diag(temp)); %#ok<SAGROW>
        pd_no_icd_index = pd_no_icd_index + 1;
    end
end

figure, imagesc(control(:,:,1)), axis equal tight, colorbar;
figure, imagesc(pd_icd(:,:,1)), axis equal tight, colorbar;
figure, imagesc(pd_no_icd(:,:,1)), axis equal tight, colorbar;

control_pd_icd = cat(3,control,pd_icd);
control_pd_no_icd = cat(3,control,pd_no_icd);


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



