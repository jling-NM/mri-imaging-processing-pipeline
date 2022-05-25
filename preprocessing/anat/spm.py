#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""@package anat

    Created on Wed Mar  3 12:08:48 2021
    
    @author: jling

    PRISMA
    
    Processing module for anat data.
    
    BIDS-dependent
 
"""

from mayerlab.preprocessing import utils
import os



def generate_tissue_classes(subj):
    """
    Invoke SPM tissue classification on T1
    :param subj: utils.Subj
    """    
    
    import subprocess, shutil
    from mayerlab.preprocessing.anat import qc
    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    if (not os.path.isfile(os.path.join(anat_path, subj.bids_file_join('spm-wm_mask.nii.gz')))):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : T1w Tissue Segmentation" )
        
        # run SPM to get segmentations
        t1_in_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
            
        t1_tmp = os.path.join(anat_path, 'tmp.T1w.nii')
        t1_tmp_gz = os.path.join(anat_path, 'tmp.T1w.nii.gz')
        
        if not os.path.isfile(t1_tmp):
            if not os.path.isfile(t1_tmp_gz):
                shutil.copyfile(t1_in_path, t1_tmp_gz)
                cmd = ("pigz -d " + t1_tmp)
                utils.Cmd(cmd).echo()
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
         
    
        if os.path.isfile(t1_tmp):
            utils.prettyOut("ANAT: generate spm tissue classes")
         
            cmd = ("matlab -nodisplay -nosplash -nojvm -singleCompThread -r \"spm_new_segment('" + t1_tmp + "')\"")
            utils.Cmd(cmd).echo()
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
        
            try:
                os.remove(t1_tmp)
            except:
                raise
            

        
        # Rename SPM segmentation output and create a binary mask
        tissue_dict = dict(gm='1',wm='2',csf='3')
        for key in tissue_dict:
            cls_file = os.path.join(anat_path, 'c'+ tissue_dict[key] +'tmp.T1w.nii')
            prob_file = os.path.join(anat_path, subj.bids_file_join('spm-' + key + '_prob.nii'))
            os.rename(cls_file, prob_file)
            utils.Cmd("pigz " + prob_file).run()
            utils.Cmd("3drefit -space ORIG -view orig " + prob_file).run()
            
            # create a binary of the mask at .60 threshold
            mask_file = os.path.join(anat_path, subj.bids_file_join('spm-' + key + '_mask.nii.gz'))
            if not os.path.isfile(mask_file):
                
                # another way to do this and include the partial volumeing would be to assign voxels to 1
                # if they are more than a third of probability and the others are less than a third
                # With the threshold method you will get overlap between tissue assignation                        
                utils.Cmd("3dcalc -a " + prob_file + ".gz -prefix " + mask_file + " -datum short -expr 'astep(a,0.60)'" ).run()
                utils.Cmd("3drefit -space ORIG -view orig " + mask_file).run()
                
                
            
        utils.purge_paths(os.path.join(anat_path, 'mtmp.T1w.nii'), os.path.join(anat_path, 'tmp.T1w_seg8.mat'))
        # take qc snapshots
        qc.snapshot_tissue_classes(subj)

                
        
                
        
        
def get_tissue_probability_file(subj, tissue='gm'):
    
    fp = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-'+tissue+'_prob.nii.gz'))
    if os.path.isfile(fp):
        return fp
    else:
        raise ValueError("No SPM probability file for " + tissue + " was found")

        
        
                
def get_tissue_mask_file(subj, tissue='gm'):
    
    fp = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-'+tissue+'_mask.nii.gz'))
    if os.path.isfile(fp):
        return fp
    else:
        raise ValueError("No SPM mask file for " + tissue + " was found")

                


def qa_tissue_classes(subj,study):

    """
    QA display of SPM tissue classification final masks.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    # afni driver doesn't appear to respond to SET_SESSION, go to directory in OS
    os.chdir(subj.derivatives_path(bidsType='anat'))
        
    # let afni know not to detach
    os.environ["AFNI_DETACH"] = "NO"

    # run the viewing command
    utils.Cmd('afni -noplugins \
    -com "OPEN_WINDOW A.coronalimage geom=320x320+889+457" \
    -com "SET_XHAIRS axialimage.OFF" \
    -com "SET_DICOM_XYZ 12 0 29" \
    -com "SWITCH_UNDERLAY '+ subj.bids_file_join('T1w.nii.gz') + '" \
    -com "SWITCH_OVERLAY '+ subj.bids_file_join('spm-csf_mask.nii.gz') + '" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 3000" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 300" \
    -com "SWITCH_OVERLAY '+ subj.bids_file_join('spm-wm_mask.nii.gz') + '" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 3000" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 300" \
    -com "SWITCH_OVERLAY '+ subj.bids_file_join('spm-gm_mask.nii.gz') + '" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +" ').run()

    
    os.environ["AFNI_DETACH"] = "YES"
    
