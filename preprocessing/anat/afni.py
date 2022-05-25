#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""@package anat

    Created on Wed Mar  3 12:08:48 2021
    
    @author: jling

    PRISMA
    
    Processing module for anat data.
    
    BIDS-dependent
 
"""

import os, shutil
from mayerlab.preprocessing import utils

def crop_t1(subj):
    """
    we want to crop the T1 before using for anything
    """

    # pull original image from rawdata source
    i_path = subj.bids_file_path(bidsType='anat', bidsLabel='T1w', extension='nii.gz')
    # output to derivaties path
    anat_path = subj.derivatives_path(bidsType='anat')
    utils.make_path(anat_path)
    o_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    
    # crop
    if not os.path.isfile(o_path):
        utils.prettyOut("ANAT : Crop T1w FOV")
        utils.Cmd("3dautobox -input {} -prefix {}".format(i_path, o_path)).run()
    
    
    
    
    
def warp_t1_to_MNI(subj):
    """
    Datasets ~2~
    
      Suppose the -prefix is 'sub007'.
      Then the outputs from this script will be"
    
      anatDO.sub007.nii       = deobliqued version of original dataset;
                                (*only if* using '-deoblique' opt);
      anatU.sub007.nii        = intensity uniform-ized original dataset
                                (or, if '-unifize_off' used, a copy of orig dset);
      anatUA.sub007.nii       = anisotropically smoothed version of the above
                                (or, if '-aniso_off' used, a copy of anatU.*.nii)
      anatUAC.sub007.nii      = ceiling-capped ver of the above (at 98%ile of 
                                non-zero values)
                                (or, if '-ceil_off' used, a copy of anatUA.*.nii)
    
      anatS.sub007.nii        = first pass skull-stripped original dataset
                                (or, if '-init_skullstr_off' used, a copy of
                                anatUAC.*.nii);
      anatSS.sub007.nii       = second pass skull-stripped original dataset;
                                * note that anatS and anatSS are 'original'
                                  in the sense that they are aligned with
                                  the input dataset - however, they have been
                                  unifized and weakly smoothed: they are
                                  stripped versions of anatUAC; if you want
                                  a skull-stripped copy of the input with
                                  no other processing, use a command like
                                    3dcalc -a INPUTDATASET        \
                                           -b anatSS.sub007.nii   \
                                           -expr 'a*step(b)'      \
                                           -prefix anatSSorig.sub007.nii
    
      anatQQ.sub007.nii       = skull-stripped dataset nonlinearly warped to
                                the base template space;
      anatQQ.sub007.aff12.1D  = affine matrix to transform original dataset
                                to base template space;
      anatQQ.sub007_WARP.nii  = incremental warp from affine transformation
                                to nonlinearly aligned dataset;
    
      * The .aff12.1D and _WARP.nii transformations need to be catenated to get
        the full warp from orginal space to the base space; example:
    
          3dNwarpApply -nwarp 'anatQQ.sub007_WARP.nii anatQQ.sub007.aff12.1D' ...
    
    QC images ~2~
    
      AMsub007.jpg            = 3x3 snapshot image of the anatQQ.sub007.nii
                                dataset with the edges from the base template
                                overlaid -- to check the alignment;
      MAsub007.jpg            = similar to the above, with the roles of the
                                template and the anatomical datasets reversed.
      QC_anatQQ.sub007.jpg    = like AM*.jpg, but 3 rows of 8 slices
      QC_anatSS.sub007.jpg    = check skullstripping in orig space: ulay is
                                input dset, and olay is mask of
                                skullstripped output (anatSS* dset)
    
      init_qc_00_overlap_uinp_obase.jpg
        o [ulay] original source dset
          [olay] original base dset
        o single image montage to check initial overlap of source and base,
          ignoring any obliquity that might be present (i.e., the way AFNI 
          GUI does by default, and also how alignment starts)
        o if initial overlap is not strong, alignment can fail or
          produce weirdness
        o *if* either dset has obliquity, then an image of both after 
          deobliquing with 3dWarp is created (*DEOB.jpg), and a text file
          about obliquity is also created (*DEOB.txt).
    
            
    """
    try:
        
        # we want to crop the T1 before using for anything
        crop_t1(subj)

        
        anat_path = subj.derivatives_path(bidsType='anat')
        i_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
        qc_path = os.path.join(subj.derivatives_path(bidsType='anat'), 'qc')
        
        if i_path:
            try:
                
                os.environ["OMP_NUM_THREADS"] = "8"
                os.environ["AFNI_ENVIRON_WARNINGS"] = "NO"
                
                # 20210304 - JL
                # while making images SSwarper using 2dcat which is looking for shared libaries
                # let it find them
                os.environ["LD_LIBRARY_PATH"] = utils.Env.afni_path
                
                
                # TODO: skulstripping now so good on single subject; if persists try:
                if not os.path.isfile(os.path.join(anat_path, subj.bids_file_join('T1w_MNI.nii.gz'))):

                    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : T1w to Anat Registration" )

                    #utils.Cmd("@SSwarper -input "+i_path+" -base MNI152_2009_template_SSW.nii.gz -subid "+subj.bids_id+" -odir "+qc_path+" -tmp_name_nice").run()
                    utils.Cmd("@SSwarper -SSopt '-push_to_edge -blur_fwhm 2' -input "+i_path+" -base MNI152_2009_template_SSW.nii.gz -subid "+subj.bids_id+" -odir "+qc_path+" -tmp_name_nice").run()

                    # compress ouput to save space
                    utils.Cmd("pigz " + qc_path + "/*.nii").run()
                    
                    # move and rename output
                    shutil.move(os.path.join(qc_path, 'anatQQ.'+subj.bids_id+'.nii.gz'), os.path.join(anat_path, subj.bids_file_join('T1w_MNI.nii.gz')) )
                    shutil.move(os.path.join(qc_path, 'anatQQ.'+subj.bids_id+'_WARP.nii.gz'), os.path.join(anat_path, subj.bids_file_join('T1w-to-MNI.WARP.nii.gz')))
                    shutil.move(os.path.join(qc_path, 'anatQQ.'+subj.bids_id+'.aff12.1D'), os.path.join(anat_path, subj.bids_file_join('T1w-to-MNI.aff12.1D')))
                    shutil.move(os.path.join(qc_path, 'anatSS.'+subj.bids_id+'.nii.gz'), os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz')))
                
            except:
                
                # clean tmp files on crash
                utils.purge_paths(qc_path+"/junk_ssw*")
                raise
                
        else:
            raise ValueError("file not found")
        
    except:
        raise
        
        
        
        
        
        
def skullstrip(subj, orig_vol=False):
    
    """
    laptop_20753 skullstrip T1 methods
    
    """    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    i_file = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    o_file = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz'))
    
    if not os.path.isfile(o_file):

        # we want to crop the T1 before using for anything
        crop_t1(subj)

        
        import uuid
        
        utils.prettyOut("SkullStrip anat...")

        tmp_dir_path = os.path.join(anat_path, str(uuid.uuid4()))
        utils.make_path(tmp_dir_path)
        os.chdir(tmp_dir_path)
        
        try:
            
            # unifize -clfrac 0.4
            utils.Cmd("3dUnifize -GM -clfrac 0.1 -Urad 30 -input " + i_file + " -prefix unif.nii.gz").run()
            
            # skullstrip
            utils.Cmd("3dSkullStrip -input unif.nii.gz -prefix sksp.tmp.nii.gz -debug 1 -ld 23 -niter 777 -shrink_fac_bot_lim 0.75 -exp_frac 0.0666 -orig_vol -pushout -use_edge -touchup").run()
            
            # smooth and fill the sksp a bit
            utils.Cmd("3dmask_tool -input sksp.tmp.nii.gz -union -dilate_input 3 -2 -prefix sksp.mask.nii.gz").run()
    
            # mask our original or unifized anat for a final sksp volume
            if orig_vol:
                utils.Cmd("3dcalc -a sksp.mask.nii.gz -b " + i_file + " -expr 'step(a)*b' -prefix " + o_file).run()
            else:
                utils.Cmd("3dcalc -a sksp.mask.nii.gz -b unif.nii.gz -expr 'step(a)*b' -prefix " + o_file).run()

        except:
            raise
            
        finally:
            # clean up
            utils.purge_paths(tmp_dir_path)
        
        
        
        
        
def affine_align_to_Template(subj):
    
    """
    laptop_20753 affine align T1 to Template (TT_N27)
    
    """        
    
    anat_path = subj.derivatives_path(bidsType='anat')
    o_file = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.AFF.nii.gz'))
    
    if not os.path.isfile(o_file):
        try:
            utils.prettyOut("affine_align_T1_to_Template anat...")
            
            #    sad limitation of this auto_tlrc
            os.chdir(anat_path)
            utils.Cmd("@auto_tlrc -no_ss -suffix 'NONE' -rmode quintic -base TT_N27+tlrc. -input " + subj.bids_file_join('T1w_SKSP.nii.gz') + " -init_xform AUTO_CENTER -prefix " + subj.bids_file_join('T1w_SKSP.TLRC.AFF.nii.gz')).run()
                
        except:
            raise
            
        finally:
            # clean up
            utils.purge_paths(os.path.join(anat_path, '*.Xat.1D'), os.path.join(anat_path, '*.nii.Xaff12.1D'), os.path.join(anat_path, '*.nii_WarpDrive.log') )

    


        
def warp_t1_to_TLRC(subj):
    
    """
    laptop_20753 Qwarp T1 to Template (TT_N27)
    
    """    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    
    affine_file = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.AFF.nii.gz'))
    aff_3x4_tx = os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC.aff12.1D'))
    final_WARP_file = os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC.NWARP.nii.gz'))
    anat_qwarp_file = os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC.nii.gz'))
    template_file = utils.Env.afni_path + "/TT_N27+tlrc"
    o_file = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.nii.gz'))


    #    final out exists?
    if not os.path.isfile(o_file):
        
        try:
            #    look for 3x4 T1 to tlrc matrix, if not there, create         
            if not os.path.isfile(final_WARP_file):
                
                if not os.path.isfile(aff_3x4_tx):
                    # get affine starting point
                    affine_align_to_Template(subj)
                    
                    utils.prettyOut("Create 3x4 T1 to tlrc transform")
                    utils.Cmd("cat_matvec " + affine_file + "::WARPDRIVE_MATVEC_FOR_000000 >  " + aff_3x4_tx).run()     
                
                
                utils.prettyOut("Create Nonlinear WARP")
                os.environ['OMP_NUM_THREADS'] = '6'
                utils.Cmd("3dQwarp -nodset -prefix " + anat_qwarp_file + " -useweight -pblur 0.09 0.09 -Qfinal -workhard:5:8 -base " + template_file + " -source " + affine_file ).run()
                # change output filename
                os.rename(os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC_WARP.nii.gz')),
                          os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC.WARP.nii.gz')))
                

                #    create final native T1 to NL tlrc WARP
                utils.prettyOut("create final native T1 to NL tlrc WARP")
                utils.Cmd("3dNwarpCat -prefix " + final_WARP_file + " -warp1 " + os.path.join(anat_path, subj.bids_file_join('T1w-to-TLRC.WARP.nii.gz')) + " -warp2 " + aff_3x4_tx ).run()
            
            #    apply final nl transform to original anat.sksp > anat.sksp.qwarp+tlrc
            utils.prettyOut("apply final nl transform to T1w")
            utils.Cmd("3dNwarpApply -interp 'wsinc5' -prefix "+o_file+" -nwarp " + final_WARP_file + " -source " + os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz')) + " -master " + template_file ).run()

        except:
            raise
            
        finally:
            utils.purge_paths(affine_file)
            utils.purge_paths(os.path.join(anat_path, 'pre.*'))
            utils.prettyOut("Completed T1 to tlrc")





def set_header_space(space='', f_path=''):
    """
    Alters file header for custom AFNI space indicator so file can be viewed in AFNI GUI
    
    "if present (sform_code > 0), is to be used to give the location of the voxels in some standard space.
    The sform_code indicates which standard space is present." - THE NIFTI-1 DATA FORMAT
    """

    if len(space) == 0:
        raise ValueError('set_afni_space needs some SPACE')
        
    if len(f_path) == 0:
        raise ValueError('set_afni_space needs a file path')
        
    if os.path.isfile(f_path):
        
        if space.upper() == 'MNI':
            sform_code = 4
        if space.upper() == 'TLRC':
            sform_code = 3
        if space.upper() == 'ORIG':
            sform_code = 1
            
        utils.prettyOut("ANAT : Set sform code to {}".format(sform_code))
        # i didn't like how AFNI was messing with precision the transform matrix so i use fslorient
        #utils.Cmd("3drefit -space {} {}".format(space, f_path)).run()
        utils.Cmd("fslorient -setsformcode {} {}".format(sform_code, f_path)).run()
        


