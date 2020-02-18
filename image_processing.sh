#!/bin/bash
subjects=`ls -d */`
dir=`pwd`
for s in $subjects
do
  code=`echo $s | cut -d"/" -f1`
  cd "${dir}/${code}"
  echo "$code holefilling.."
  thresh_it=`fslstats T1/T1_WM_mask.nii.gz -P 15`
  fslmaths T1/T1_WM_mask.nii.gz -thr $thresh_it -bin T1/T1_WM_mask_bin
  fslmaths T1/T1_WM_mask_bin.nii.gz -kernel 3D -dilM -ero T1/T1_WM_mask_bin_holefilled.nii.gz
  fslmaths T1/T1_brain_mask.nii.gz -kernel 3D -dilM -ero T1/T1_brain_mask_holefilled.nii.gz
  fslmaths T1/T1.nii.gz -mul T1/T1_brain_mask_holefilled.nii.gz T1/T1_brainMasked

  echo "Label conversion.."
  labelconvert T1/GM_parcellation.nii.gz ../FreeSurferColorLUT.txt ../fs_default_Bstem.txt T1/GM_parcellation_85labels.nii.gz

  echo "Registrering T1 on B0.."
  mkdir -p 'Registrations/Reg_T1onB0'
  fslroi DTI/data_preproc.nii.gz DTI/data_preproc_b0.nii.gz 0 1
  flirt -in DTI/data_preproc_b0.nii.gz -ref T1/T1_brainMasked.nii.gz -dof 6 -omat Registrations/Reg_T1onB0/temp.mat
  flirt -in DTI/data_preproc_b0.nii.gz -ref T1/T1_brainMasked.nii.gz -dof 6 -cost bbr \
    -wmseg T1/T1_WM_mask_bin_holefilled.nii.gz -init Registrations/Reg_T1onB0/temp.mat \
    -omat Registrations/Reg_T1onB0/b0_T1.mat -out Registrations/Reg_T1onB0/b0_T1
  convert_xfm -omat Registrations/Reg_T1onB0/T1_b0.mat -inverse Registrations/Reg_T1onB0/b0_T1.mat
  flirt -in T1/T1_brainMasked.nii.gz -ref DTI/data_preproc_b0.nii.gz -applyxfm -init Registrations/Reg_T1onB0/T1_b0.mat \
    -dof 6 -out Registrations/Reg_T1onB0/T1onB0

  echo "Registrering WMM,BM,GMP on B0.."
  mkdir 'Registrations/Reg_WMonB0'
  mkdir 'Registrations/Reg_brainMaskonB0'
  mkdir 'Registrations/Reg_GMonB0'
  flirt -in T1/T1_WM_mask_bin_holefilled.nii.gz -ref DTI/data_preproc_b0.nii.gz -applyxfm -init Registrations/Reg_T1onB0/T1_b0.mat \
    -interp nearestneighbour -dof 6 -out Registrations/Reg_WMonB0/WMonB0
  flirt -in T1/T1_brain_mask_holefilled.nii.gz -ref DTI/data_preproc_b0.nii.gz -applyxfm -init Registrations/Reg_T1onB0/T1_b0.mat \
    -interp nearestneighbour -dof 6 -out Registrations/Reg_brainMaskonB0/brainMaskonB0
  flirt -in T1/GM_parcellation_85labels.nii.gz -ref DTI/data_preproc_b0.nii.gz -applyxfm -init Registrations/Reg_T1onB0/T1_b0.mat \
    -interp nearestneighbour -dof 6 -out Registrations/Reg_GMonB0/GMonB0
  # why -interp nearestneighbour? The main difference is for GM_parcellation_85labels, used in connectome generation. Trilinear = blurred matrix, nearestneighbour = nice isolated sections

  echo "Computing CSD and peaks selection.."
  mkdir -p 'Tractography/CSD'
  dwi2response tournier DTI/data_preproc.nii.gz Tractography/CSD/response.txt -fslgrad DTI/data_preproc.bvecs DTI/data_preproc.bvals \
    -mask Registrations/Reg_brainMaskonB0/brainMaskonB0.nii.gz -lmax 6
  dwi2fod csd DTI/data_preproc.nii.gz Tractography/CSD/response.txt Tractography/CSD/FODsh.mif -fslgrad DTI/data_preproc.bvecs DTI/data_preproc.bvals \
    -lmax 6 -mask Registrations/Reg_brainMaskonB0/brainMaskonB0.nii.gz
  sh2peaks -threshold 0.2 Tractography/CSD/FODsh.mif Tractography/CSD/FODsh_peaks_0.2.mif

  echo "Generation streamlines.."
  tckgen Tractography/CSD/FODsh_peaks_0.2.mif Tractography/CSD/fibers_det_1000000.tck -algorithm FACT -seed_image Registrations/Reg_WMonB0/WMonB0.nii.gz \
    -mask Registrations/Reg_WMonB0/WMonB0.nii.gz -select 1000000 -force

  echo "Computing connectome.."
  mkdir 'Connectome'
  tck2connectome Tractography/CSD/fibers_det_1000000.tck Registrations/Reg_GMonB0/GMonB0.nii.gz Connectome/${code}_connectome.csv \
    -nthreads 4 -force #assignment_radial_search: default
done

# If we use trilinear GM interpolation, when we visualize the connectome as a graph we can see that there is something wrong
# mrconvert -datatype uint32 Registrations/Reg_GMonB0/GMonB0.nii.gz parcInt.nii.gz
# mrview Registrations/Reg_T1onB0/T1onB0.nii.gz -connectome.init parcInt.nii.gz -connectome.load Connectome/HC_3104_connectome.csv
