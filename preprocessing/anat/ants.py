#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 11:37:16 2021

@author: jling
"""


import os
from mayerlab.preprocessing import utils
#import ants  # this loads really slow at this time so avoiding py interface
    
    

def brain_extraction(subj):
    """
    What it says...
    """
    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    i_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    o_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    
    if not os.path.isfile(o_path):
        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "22"
        
        utils.prettyOut("ANAT : Brain Extraction")
        utils.Cmd('antsBrainExtraction_fmriprep.sh \
        -d 3 \
        -a {} \
        -e /export/research/analysis/human/amayer/shared/apps/mindboggle/OASIS-30_Atropos_template/T_template0.nii.gz \
        -m /export/research/analysis/human/amayer/shared/apps/mindboggle/OASIS-30_Atropos_template/T_template0_BrainCerebellumProbabilityMask.nii.gz \
        -f /export/research/analysis/human/amayer/shared/apps/mindboggle/OASIS-30_Atropos_template/T_template0_BrainCerebellumRegistrationMask.nii.gz \
        -o {}'.format( i_path, os.path.join(anat_path, subj.bids_file_join('T1w_')) )).run()
        
    
        # rename transform output for consistency
        os.rename(os.path.join(anat_path, subj.bids_file_join('T1w_BrainExtractionPrior0GenericAffine.mat')), os.path.join(anat_path, subj.bids_file_join('T1w-to-mni_OASIS-30.aff.transform.mat')))
        # rename masked data
        os.rename(os.path.join(anat_path, subj.bids_file_join('T1w_BrainExtractionBrain.nii.gz')), os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz')))
        # rename brain mask
        os.rename(os.path.join(anat_path, subj.bids_file_join('T1w_BrainExtractionMask.nii.gz')), os.path.join(anat_path, subj.bids_file_join('T1w_brain_mask.nii.gz')))
        # this removes other leftovers which are mostly the ants tissue segmentation which is problematic. We use SPM
        utils.purge_paths(os.path.join(anat_path, subj.bids_file_join('T1w_BrainExtraction*')))




def denoise(subj, i_path="", o_path=""):
    """
    A denoising method for MRI
    """
    
    anat_path = subj.derivatives_path(bidsType='anat')
    if len(i_path) == 0:
        i_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    if len(o_path) == 0:
        o_path = os.path.join(anat_path, subj.bids_file_join('T1w-den.nii.gz'))
    
    
    if not os.path.isfile(o_path):
        utils.prettyOut("ANAT : Denoise T1w")
        utils.Cmd('DenoiseImage \
                    --image-dimensionality 3 \
                    --input-image {} \
                    --noise-model Rician \
                    --shrink-factor 1 \
                    --patch-radius 1 \
                    --search-radius 2 \
                    --output {} \
                    --verbose 1'.format(i_path, o_path)).run()
                    


    

def affine_initializer(subj, template_key='mni_icbm152_nlin_asym_2009c'):
    anat_path = subj.derivatives_path(bidsType='anat')
    i_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    o_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'.aff.transform.mat'))
    
    if not os.path.isfile(o_path):
        utils.prettyOut("ANAT : Affine Registration")
        utils.Cmd('antsAffineInitializer 3 {} {} {} 15.000000 0.100000 0 10'.format(utils.Env.template_lib[template_key]['brain'], i_path, o_path)).run()
        




    
def warp_t1_to_MNI(subj, template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Registration step to create composite transform (*.h5)
    """
    
    anat_path = subj.derivatives_path(bidsType='anat')
    i_path_tx = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'.aff.transform.mat'))
    i_path_brain = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    o_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_'))
    o_path_tx = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))
    
    if not os.path.isfile(o_path_tx):
        
        # affine init prior to NL reg
        affine_initializer(subj, template_key=template_key)

        
        utils.prettyOut("ANAT : Nonlinear Registration")
        
        os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "22"
        
        utils.Cmd('antsRegistration \
                    --collapse-output-transforms 1 \
                    --dimensionality 3 \
                    --float 0 \
                    --initial-moving-transform [ {0}, 0 ] \
                    --initialize-transforms-per-stage 0 \
                    --interpolation LanczosWindowedSinc \
                    --output [ {1} ] \
                    --transform Rigid[ 0.05 ] \
                    --metric Mattes[ {2}, {3}, 1, 56, Regular, 0.25 ] \
                    --convergence [ 100x100, 1e-06, 20 ] \
                    --smoothing-sigmas 2.0x1.0vox \
                    --shrink-factors 2x1 \
                    --use-estimate-learning-rate-once 1 \
                    --use-histogram-matching 1 \
                    --transform Affine[ 0.08 ] \
                    --metric Mattes[ {2}, {3}, 1, 56, Regular, 0.25 ] \
                    --convergence [ 100x100, 1e-06, 20 ] \
                    --smoothing-sigmas 1.0x0.0vox \
                    --shrink-factors 2x1 \
                    --use-estimate-learning-rate-once 1 \
                    --use-histogram-matching 1 \
                    --transform SyN[ 0.1, 3.0, 0.0 ] \
                    --metric CC[ {2}, {3}, 1, 4, None, 1 ] \
                    --convergence [ 100x70x50x20, 1e-06, 10 ] \
                    --smoothing-sigmas 3.0x2.0x1.0x0.0vox \
                    --shrink-factors 8x4x2x1 \
                    --use-estimate-learning-rate-once 1 \
                    --use-histogram-matching 1 \
                    --winsorize-image-intensities [ 0.005, 0.995 ] \
                    --write-composite-transform 1 \
                    --verbose'.format( i_path_tx, o_path, utils.Env.template_lib[template_key]['brain'], i_path_brain )).run()
        
        
        utils.purge_paths(os.path.join(anat_path, subj.bids_file_join('T1w_N4*.nii.gz')))


        
        
        
def apply_spatial_normalization(subj, template_key='mni_icbm152_nlin_asym_2009c', i_path="", o_path=""):

    tx_path = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))

    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "22"
    
    utils.Cmd('antsApplyTransforms \
        --default-value 0 \
        --dimensionality 3 \
        --float 0 \
        --input {} \
        --interpolation LanczosWindowedSinc \
        --output {} \
        -r {} \
        -t {}'.format(i_path, o_path, utils.Env.template_lib[template_key]['brain'], tx_path)).run()
    
    
    # update NIFTI header
    from mayerlab.preprocessing.anat import afni
    afni.set_header_space(space='MNI', f_path=o_path)


        
