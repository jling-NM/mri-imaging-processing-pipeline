# -*- coding: utf-8 -*-

"""@package asl
Preprocessing module for arterial spin label imaging

Created on Wed Dec  2 13:28:49 2015

@author: jling

"""

import os
from mayerlab.preprocessing import utils

class pCASLMbMd:
    
    class pipelines:
            
        def standard(subj):
            """
            #=====================================================================================
            # Xiufeng Li, Ph.D.
            # CMRR, UoM
            # Minneapolis, MN, USA
            # 2015-12-10
            #----------------------------------
            
            #============================================================================
            # Multi-Delay PCASL Configuration Tables
            # 1. Three tables: DELAYS, MEASUREMENTS and optional INTERLEAVED ACQUISITION
            # 
            # 2. Max. # of table elements = xxx_Element0 = 10!
            # 
            # 3. Minimal MEASUREMENTS value > 2 for each element (delay)!
            #============================================================================
            
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #-Table 1: DELAY TIMES (us) in order;
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            [MD_ASL_DELAYS]
            MD_ASL_DELAYS_Element0   = 10
            
            #------------------------------------------------------
            #[High Resolution]
            #------------------------------------------------------
            MD_ASL_DELAYS_Element1  =  200000
            MD_ASL_DELAYS_Element2  =  700000
            MD_ASL_DELAYS_Element3  = 1200000
            MD_ASL_DELAYS_Element4  = 1700000
            MD_ASL_DELAYS_Element5  = 2200000
            
            MD_ASL_DELAYS_Element6  = 2200000
            MD_ASL_DELAYS_Element7  = 2200000
            MD_ASL_DELAYS_Element8  = 2200000
            MD_ASL_DELAYS_Element9  = 2200000
            MD_ASL_DELAYS_Element10 = 2200000
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #-Table 2: MEASUREMENTS NUMBER in order
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            [MD_ASL_MEASUREMENTS]
            MD_ASL_MEASUREMENTS_Element0  = 10
            
            #------------------------------------------------------
            #[High Resolution]:  NEX = 86 + 2 Noise + 2 M0 = UI 90 for 5 mins 30 s!
            #Proposed for HCP lifespan study
            #------------------------------------------------------
            MD_ASL_MEASUREMENTS_Element1 = 12	
            MD_ASL_MEASUREMENTS_Element2 = 12	
            MD_ASL_MEASUREMENTS_Element3 = 12	
            MD_ASL_MEASUREMENTS_Element4 = 20
            MD_ASL_MEASUREMENTS_Element5 = 30
            
            MD_ASL_MEASUREMENTS_Element6 = 30
            MD_ASL_MEASUREMENTS_Element7 = 40
            MD_ASL_MEASUREMENTS_Element8 = 50
            MD_ASL_MEASUREMENTS_Element9 = 60
            MD_ASL_MEASUREMENTS_Element10= 70
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #-Table 3: Interleaved Ascending and Decending Acquisition Flag:
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            [MD_ASL_INTERLEAVED_ORDER]
            MD_ASL_INTERLEAVED_ORDER_Element0  = 10
            
            MD_ASL_INTERLEAVED_ORDER_Element1  = 0
            MD_ASL_INTERLEAVED_ORDER_Element2  = 0         
            MD_ASL_INTERLEAVED_ORDER_Element3  = 0
            MD_ASL_INTERLEAVED_ORDER_Element4  = 0
            MD_ASL_INTERLEAVED_ORDER_Element5  = 0
            
            MD_ASL_INTERLEAVED_ORDER_Element6  = 0
            MD_ASL_INTERLEAVED_ORDER_Element7  = 0
            MD_ASL_INTERLEAVED_ORDER_Element8  = 0
            MD_ASL_INTERLEAVED_ORDER_Element9  = 0
            MD_ASL_INTERLEAVED_ORDER_Element10 = 0
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<        
    
    
    
            # convert PLDs(ms) to TIs(sec)
            [TI = labeling duration + post labelling delay.]
            1500 + 200 = 1700 ms ( 1.7 s)
            1500 + 700 = 2200 ms ( 2.2 s)
            1500 + 1200 = 2700 ms ( 2.7 s)
            1500 + 1700 = 3200 ms ( 3.2 s)
            1500 + 2200 = 3700 ms ( 3.7 s)
            
            1.7,2.2,2.7,3.2,3.7
            """
            
            
            perf_raw_path = subj.bids_modality_path(bidsType='perf')
            if not os.path.isdir(perf_raw_path):
                raise ValueError("Missing input:\n"+perf_raw_path)
        
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(perf_o_path)
            
                    
            ### Combine fmaps into single file
            topup_i_file = os.path.join(perf_o_path, subj.bids_file_join("perf-per.nii.gz"))
            fmap_path = subj.bids_modality_path(bidsType='fmap')
            ap_file = os.path.join(fmap_path, subj.bids_file_join("acq-pcaslmbmd_dir-ap_epi.nii.gz"))
            pa_file = os.path.join(fmap_path, subj.bids_file_join("acq-pcaslmbmd_dir-pa_epi.nii.gz"))
        
            if not os.path.isfile(topup_i_file):    
                if os.path.isfile(ap_file) and os.path.isfile(pa_file):
                    utils.Cmd("fslmerge -t {} {} {} ".format(topup_i_file, ap_file, pa_file) ).run()            
            
            
        
            ### Estimate distortion from fmaps
            dewarped_file = os.path.join(perf_o_path, "topup.dewarped_fmap.nii.gz")
            if not os.path.isfile(dewarped_file):
                utils.Cmd( "topup --verbose"
                        + " --imain=" + topup_i_file
                        + " --datain=" + subj.study.paramsPath + '/asl_params.txt'
                        + " --config=b02b0.cnf"                        
                        + " --out=" + os.path.join(perf_o_path, "topup_results")
                        + " --iout=" + dewarped_file
                        + " --logout=" + os.path.join(perf_o_path, "topup_log") ).run()
            
            
    
                # save qc images
                pCASLMbMd.qc.snapshot_topup(subj, topup_i_file, dewarped_file)

                
                
            
            ### Create a mask from mean dewarped fmaps
            dewarped_Tmean_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean.nii.gz")
            if not os.path.isfile(dewarped_Tmean_file): 
                utils.Cmd("fslmaths {} -Tmean {}".format(dewarped_file, dewarped_Tmean_file)).run()
            
            dewarped_Tmean_mask_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean")
            if not os.path.isfile(dewarped_Tmean_mask_file + '_mask.nii.gz'): 
                utils.Cmd("bet {} {} -f 0.45 -g 0.10 -n -m".format(dewarped_Tmean_file, dewarped_Tmean_mask_file)).run()
                
                
                
            ### Spatial normalization from mean dewarped fmaps
            pCASLMbMd.spatial_normalization(subj)
            
            
            asl_i_file = os.path.join(perf_raw_path, subj.bids_file_join("acq-pcaslmbmd_asl.nii.gz"))
            asl_dc_file = os.path.join(perf_o_path, subj.bids_file_join("acq-pcaslmbmd_asl_dc.nii.gz"))
            m0_file = os.path.join(perf_o_path, subj.bids_file_join("acq-pcaslmbmd_M0.nii.gz"))
            asl_file = os.path.join(perf_o_path, subj.bids_file_join("acq-pcaslmbmd_asl-meandiff.nii.gz"))
            
            
            ### Apply distortion to raw asl data
            if (not os.path.isfile(m0_file)) or (not os.path.isfile(asl_file)):
                if not os.path.isfile(asl_dc_file):
                    utils.Cmd("applytopup " + 
                              " --verbose " +
                              " --imain=" + asl_i_file +
                              " --inindex=1" + 
                              " --topup=" + os.path.join(perf_o_path, 'topup_results_fieldcoef') +
                              " --datain=" + subj.study.paramsPath  + "/asl_params.txt"
                              " --method=jac " +
                              " --interp=spline" + 
                              " --out=" + asl_dc_file).run()
            


            ### Extract M0 images to mean M0
            if not os.path.isfile(m0_file):

                # extract/register M0 to dewarped_Tmean_file
                utils.Cmd("3dvolreg -overwrite -base {0} -prefix {2} {1}'[88,89]'".format(dewarped_Tmean_file, asl_dc_file, os.path.join(perf_o_path, "tmp_M0.nii.gz"))).run()
                # get mean M0
                utils.Cmd("fslmaths {} -Tmean {} -odt short".format(os.path.join(perf_o_path, "tmp_M0.nii.gz"), m0_file)).run()
                # save qc images
                pCASLMbMd.qc.snapshot_calibration_input(subj, os.path.join(perf_o_path, "tmp_M0.nii.gz"), m0_file)
                # clean up
                utils.purge_paths(os.path.join(perf_o_path, "tmp_M0.nii.gz"))
            
            
            
            
            ### Create ASL pairs mean diff file
            if not os.path.isfile(asl_file):            
                # Extract ASL pairs and register to dewarped fmap
                utils.Cmd("3dvolreg -overwrite -base {} -prefix {} {}'[0..85]'".format(dewarped_Tmean_file, os.path.join(perf_o_path, "tmp_pairs.reg.nii.gz"), asl_dc_file)).run()
                
                # get mean diff WRT TIs and blocks
                utils.Cmd("asl_file --data={} --ntis=5 --rpts=12,12,12,20,30 --ibf=tis --obf=tis --iaf=tc --diff --mean={}".format(
                    os.path.join(perf_o_path, "tmp_pairs.reg.nii.gz"),
                    asl_file
                    )).run()
                
                # clean up
                utils.purge_paths(os.path.join(perf_o_path, "tmp_pairs*.nii.gz"))                
                
                
            # remove asl file
            utils.purge_paths(asl_dc_file)
            
            
            
            
            ### get gm and wm tissue probability into native space for CBF quantification
            from mayerlab.preprocessing.anat import spm
            gm_file = spm.get_tissue_probability_file(subj, tissue='gm')
            wm_file = spm.get_tissue_probability_file(subj, tissue='wm')
            gm_asl_file = os.path.join(perf_o_path, 'gm-prob.nii.gz')
            wm_asl_file = os.path.join(perf_o_path, 'wm-prob.nii.gz')
            dewarped_Tmean_mask_file = os.path.join(perf_o_path, 'topup.dewarped_fmap_Tmean_mask.nii.gz')
            
            
            ### Estimate CBF
            cbf_path = os.path.join(perf_o_path, 'cbf')
            if not os.path.isfile(os.path.join(cbf_path, 'native_space', 'perfusion_calib.nii.gz')):

                # quantify CBF with spatial regularization enabled
                utils.Cmd("oxford_asl -i {} -o {} -m {} --spatial=on --cmethod voxel -c {} --tis 1.7,2.2,2.7,3.2,3.7 --casl --bolus 1.5 --bat 1.3 --t1 1.3 --t1b 1.65 --slicedt 0.06000 --sliceband 10 --alpha 0.85 --tr 8.0 ".format(
                    asl_file, cbf_path, dewarped_Tmean_mask_file, m0_file)).run()


            cbf_path = os.path.join(perf_o_path, 'cbf_pvcorr')
            if not os.path.isfile(os.path.join(cbf_path, 'native_space', 'perfusion_calib.nii.gz')):

                try:
                    pCASLMbMd.apply_t1_to_native(subj, gm_file, dewarped_Tmean_file, gm_asl_file)
                    pCASLMbMd.apply_t1_to_native(subj, wm_file, dewarped_Tmean_file, wm_asl_file)
                    
                    # quantify CBF with spatial regularization disabled and partial volume correction enabled
                    utils.Cmd("oxford_asl -i {} -o {} -m {} --spatial=off --pvcorr --pvgm {} --pvwm {} --cmethod voxel -c {} --tis 1.7,2.2,2.7,3.2,3.7 --casl --bolus 1.5 --bat 1.3 --t1 1.3 --t1b 1.65 --slicedt 0.06000 --sliceband 10 --alpha 0.85 --tr 8.0 ".format(
                        asl_file, cbf_path, dewarped_Tmean_mask_file, gm_asl_file, wm_asl_file, m0_file)).run()    
                    
                except:
                    raise
                    
                finally:
                    utils.purge_paths(gm_asl_file, wm_asl_file)




                

            
            
            # ### qc images
            pCASLMbMd.qc.snapshot_spatial_normalization(subj)
            pCASLMbMd.qc.snapshot_mask(subj)
            pCASLMbMd.qc.snapshot_cbf(subj)
            pCASLMbMd.qc.cbf_plots(subj, space='MNI')
            # pdf of images
            pCASLMbMd.qc.aggregate_qc_pdf(subj)
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF Processsing Complete")
            





    
    
    
    def spatial_normalization(subj, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
        """
        Boundary-Based Registration works better than ANTs to MNI template
        
        """
        
        anat_path = subj.derivatives_path(bidsType='anat')
        t1_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
        t1_brain_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
        t1_wmseg_path = os.path.join(anat_path, subj.bids_file_join('spm-wm_prob.nii.gz'))
        
        perf_o_path = subj.derivatives_path(bidsType='perf')
        tx_for_pref = os.path.join(perf_o_path, subj.bids_file_join('perf-to-T1.aff'))
        tx_inv_pref = os.path.join(perf_o_path, subj.bids_file_join('perf-from-T1.aff'))
        
    
        try:
            
            if not os.path.isfile(tx_for_pref + '.mat'):   
                
                try:
                    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF : Spatial Normalization" )
                                
                    # check for inputs
                    fmap_input_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean.masked.nii.gz")
                    utils.mask_data(os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean.nii.gz"),
                                    os.path.join(perf_o_path, 'topup.dewarped_fmap_Tmean_mask.nii.gz'),
                                    fmap_input_file)
                    
                    # do a standard flirt pre-alignment
                    print("   +++ initial transform")
                    utils.Cmd("flirt -ref {} -in {} -dof 12 -omat {}".format(t1_brain_path, fmap_input_file, os.path.join(perf_o_path, subj.bids_file_join('perf', 'init.mat')))).run()
            
                    # now run the bbr
                    print("   +++ BBR transform")
                    utils.Cmd("flirt -ref {} -in {} -dof 12 -cost bbr -wmseg {} -init {} -omat {} -schedule {}/etc/flirtsch/bbr.sch".format(
                        t1_brain_path, fmap_input_file, t1_wmseg_path,
                        os.path.join(perf_o_path, subj.bids_file_join('perf', 'init.mat')),
                        tx_for_pref + '.mat',
                        utils.Env.fsl_path
                        )).run()
                
                except:
                    raise
                    
                finally:
                    # remove init mat
                    utils.purge_paths(os.path.join(perf_o_path, subj.bids_file_join('perf', 'init.mat')))
            
        
        
            # also save out inverse transform
            if not os.path.isfile(tx_inv_pref + '.mat'):
                print("   +++ inverse transform")
                utils.Cmd("convert_xfm -omat {} -inverse {}".format(tx_inv_pref + '.mat', tx_for_pref + '.mat')).run()
                
            
            # convert FSL transform to ITK format for ANTs
            if not os.path.isfile(tx_for_pref + '.mat.itk.txt'):
                print("   +++ convert transform for ANTs")
                utils.Cmd(utils.Env.c3d_path + "/c3d_affine_tool -ref {} -src {} {} -fsl2ras -oitk {}".format(t1_path, fmap_input_file, tx_for_pref + '.mat', tx_for_pref + '.mat.itk.txt' )).run()
              
            
            # collapse the transformations to a final displacement field
            composite_warp_path = os.path.join(perf_o_path, subj.bids_file_join('perf-to-'+template_key+'_Composite.warp.nii.gz'))
            composite_in_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))
            if not os.path.isfile(composite_warp_path):
                print("   +++ composite transform")
                utils.Cmd("antsApplyTransforms -d 3 -o '[{},1]' -t {} -t {} -r {}".format(composite_warp_path, composite_in_path, tx_for_pref + '.mat.itk.txt', utils.Env.template_lib[template_key]['brain'])).run()
    
        
        except:
            raise
            
        finally:
            utils.purge_paths(tx_for_pref + '_fast_*', tx_for_pref + '_init*',tx_for_pref + '.nii.gz')
        
        
        
        
            
            
 
        
 
    
 
    def apply_spatial_normalization(subj, i_file, o_path=None, master_path=None, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
        """
        
        Parameters
        ----------
        subj : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        """
    
        try:
    
            if not os.path.isfile(o_path):
                
                perf_o_path = subj.derivatives_path(bidsType='perf')
                
                if not os.path.isfile(i_file):
                    raise ValueError("Input file does not exist")
                    
                # split input into path and file
                (i_path, i_file_name) = os.path.split(i_file) 
                
                # use i_path as o_path if user did not specify
                if o_path == None:
                    o_path = i_path
                    utils.make_path(o_path)
                     
                if master_path == None:
                    master_path = utils.Env.template_lib[template_key]['brain']
    
        
                utils.Cmd("antsApplyTransforms -d 3 -o {} -t {} -r {} -i {} --output-data-type short --interpolation LanczosWindowedSinc --input-image-type 0".format(
                    o_path,
                    os.path.join(perf_o_path, subj.bids_file_join('perf-to-'+template_key+'_Composite.warp.nii.gz')),
                    master_path,
                    os.path.join(i_path, i_file_name)
                    )).run()
                
                
                # update NIFTI header
                from mayerlab.preprocessing.anat import afni
                afni.set_header_space(space='MNI', f_path=os.path.join(o_path))
                
        except:
            raise
    
    
    





    
    

            
            
            
    def apply_t1_to_native(subj, t1_space_i_path, asl_space_reference_i_path, asl_space_o_path):
        """
        
        Parameters
        ----------
        subj : TYPE
            Moves from T1 to ASL native space
    
        Returns
        -------
        None.
    
        """
        try:
    
            if not os.path.isfile(asl_space_o_path):
                
                perf_path = subj.derivatives_path(bidsType='perf')
                tx_inv_pref = os.path.join(perf_path, subj.bids_file_join('perf-from-T1.aff.mat'))
                
                if not os.path.isfile(t1_space_i_path):
                    raise ValueError("Input file does not exist")
                if not os.path.isfile(tx_inv_pref):
                    raise ValueError("mat file does not exist")
                    
                utils.Cmd("flirt -in {} -ref {} -out {} -init {} -applyxfm".format(
                    t1_space_i_path,
                    asl_space_reference_i_path,
                    asl_space_o_path,
                    tx_inv_pref
                    )).run()
            
        except:
            raise
            
            
            
            
            
    def apply_native_to_t1(subj, asl_space_i_path, t1_space_reference_i_path, t1_space_o_path):
        """
        
        Parameters
        ----------
        subj : TYPE
            Moves from ASL native space to T1
    
        Returns
        -------
        None.
    
        """
        try:
    
            if not os.path.isfile(t1_space_o_path):
                
                perf_path = subj.derivatives_path(bidsType='perf')
                tx_pref = os.path.join(perf_path, subj.bids_file_join('perf-to-T1.aff.mat'))
                
                if not os.path.isfile(asl_space_i_path):
                    raise ValueError("Input file does not exist")
                if not os.path.isfile(tx_pref):
                    raise ValueError("mat file does not exist")
                    
                utils.Cmd("flirt -in {} -ref {} -out {} -init {} -applyxfm".format(
                    asl_space_i_path,
                    t1_space_reference_i_path,
                    t1_space_o_path,
                    tx_pref
                    )).run()
            
        except:
            raise


            
            
            
        
    class qc:
        
        ### snapshot of two M0 images side x side and mean M0 image to verify calibration image
        # afni scripts bulky and buggy, probably should use /export/research/analysis/human/amayer/shared/apps/fsl/bin/slicer
        def snapshot_calibration_input(subj, i_file, o_file):
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(os.path.join(perf_o_path, 'qc'))

            #utils.Cmd("@djunct_4d_imager -inset {} -prefix {}".format(i_file, 'M0_input') ).run()
            utils.Cmd("@djunct_4d_imager -inset {} -prefix {}".format(i_file, os.path.join(perf_o_path, 'qc', 'M0_input')) ).run()
            
            utils.Cmd("@chauffeur_afni \
                        -ulay "+o_file+" \
                        -prefix "+ os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', 'M0_output_qc_onescl') +" \
                        -set_xhairs OFF \
                        -montx 1 \
                        -monty 1").run()
        
        
        
        
        ### snapshot of fmaps and topup results to verify distortion correction input
        def snapshot_topup(subj, i_file, dewarped_file):
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(os.path.join(perf_o_path, 'qc'))
            
            utils.Cmd("@djunct_4d_imager -inset {} -prefix {}".format(i_file, os.path.join(perf_o_path, 'qc', 'topup_input'))).run()
            utils.Cmd("@djunct_4d_imager -inset {} -prefix {}".format(dewarped_file, os.path.join(perf_o_path, 'qc', 'topup_output'))).run()

        
        
        
        
        def snapshot_mask(subj):
            """
            create mask images into qc directory
            """
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(os.path.join(perf_o_path, 'qc'))
            os.chdir( os.path.join(perf_o_path, 'qc') )
            
            ulay_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean.nii.gz")            
            olay_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean_mask.nii.gz")
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+olay_file+" \
                    -prefix mask \
                    -set_xhairs OFF \
                    -montx 4 \
                    -monty 1 \
                    -delta_slices 4 4 2").run()
                    
                    
                    
        
        def snapshot_cbf(subj):  
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(os.path.join(perf_o_path, 'qc'))
            
            ulay_file = os.path.join(perf_o_path, "topup.dewarped_fmap_Tmean.nii.gz")
            olay_file = os.path.join(perf_o_path, 'cbf', 'native_space', 'perfusion_calib.nii.gz')
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ASL : QC CBF Snapshots" )
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+olay_file+" \
                    -prefix "+ os.path.join(perf_o_path, 'qc', subj.bids_file_join('CBF')) +" \
                    -set_xhairs OFF \
                    -montx 4 \
                    -monty 3 \
                    -delta_slices 3 3 3 \
                    -opacity 9 \
                    -cbar_topval 100 \
                    -cbar_ncolors 9 \
                    -func_range 80 -thr_olay 0.1").run()
                    
                    
                    
                    
        ### plots of CBF
        def cbf_plots(subj, space='TLRC'):
            
            """
            ALSOP 2015
            https://doi.org/10.1002/mrm.25197
        
            For evaluating the quality of ASL MRI images in clinical practice we advise the following
            checks:
            
            1.
            For PCASL scans, look for areas of low labeling efficiency. First, identify which
            arteries should have been labeled. Typically, this will include internal and external
            carotid arteries, and vertebral arteries. If an angiogram is available, this can be used
            to verify the list of labeled arteries. Checking the Circle of Willis anatomy may also
            be of use in matching vascular territories to labeled arteries. When the labeling
            efficiency is low in an artery, the entire flow territory of that artery will
            demonstrate a low calculated CBF. When a low CBF area is seen that matches an
            entire vascular territory, with no apparent compensation from other arteries, a
            labeling failure should be considered, though this does not preclude the possibility
            of truly low CBF conditions, or abnormally long ATT. Labeling failures can be
            caused by tortuous vessels or resonance offsets in the labeling plane. The former
            may be addressed by adjusting the location of the labeling plane, in which case an
            additional angiogram would be helpful. The latter is commonly caused by dental
            work, and may be suggested by signal dropouts around the teeth in other images
            from this patient. Methods to address resonance offset related labeling problems in
            PCASL are discussed above.
            
            2.
            Note the overall gray matter CBF value. Absolute CBF values obtained in gray
            matter can vary significantly, even among healthy young adults, due to natural
            inter-subject and intra-subject variations. In addition, average numbers are sensitive
            to partial volume effects and the methods used for isolating gray matter signal. As a
            general rule, gray matter CBF values from 40–100 ml/min/100ml can be normal.
            When the overall gray matter CBF value is inconsistent with the expected values
            for the patient population, consider the possibility that there is a global reduction in
            labeling efficiency, or that the PD scan used for normalization was incorrectly
            acquired or scaled. Clear contrast between gray and white matter should be present,
            and if not, may signify either poor labeling or motion artifacts.
            
            3.
            Check for motion artifacts. As a subtractive technique, ASL is motion sensitive,
            though this sensitivity is mitigated by BS as discussed above. The presence of
            signal outside of the brain, frequently recognizable as signal from layers of skin or
            fat is a clear indication of significant subject motion. When possible, it may be
            useful to check individual label/control difference images before averaging to see
            whether artifacts arise from only a minority of these difference images. If so, these
            images can be excluded from the CBF calculation. In addition, motion correction
            by means of automated image registration algorithms can be performed, though
            these may not be effective when BS is very efficient, or when applied to the label/
            control difference images, as the individual image SNR in these cases is low. When
            BS is not used or is incomplete, image registration is likely to be more effective,
            but BS is nevertheless recommended as a primary means of reducing physiological
            noise and motion artifacts. In the ideal case, prospective motion correction methods
            can be used when available to reduce motion artifacts during acquisition (73), and
            some of these methods are compatible with background suppression.
            
            4.
            Look for intravascular artifacts. Hyperintense spots and serpiginous regions often
            represent intravascular signal. When observed, it is advisable to verify that the PLD
            was appropriate for the patient (see Table 1), as a low PLD will naturally generate
            ASL signals in larger arteries. Intra-arterial signal with a correct PLD suggests that
            delivery of labeled blood to tissue is delayed, through slow flow and/or circuitous
            or collateral routes of circulation. Intra-venous ASL signal suggests that an
            arteriovenous shunt is present. Note that CBF calculations over whole brain or
            large regions of interest may still be valid in the presence of intravascular artifact as
            long as flow crushing gradients were not used.
            
            5. 
            Check the borderzone (watershed) regions. The borderzone or watershed areas are
            at the more distal portions of each vascular territory, and will naturally have a
            longer ATT than other portions of the territory. Note that it is possible for low ASL
            signal in these regions to represent long ATT rather than low CBF, and an
            additional scan with longer PLD may help to distinguish between these two
            possibilities. An example of this effect is shown in Figure 9.
            """
        
            import nibabel as nib
            import os
            import seaborn as sns
            import matplotlib.pyplot as plt
        
            utils.prettyOut(subj.ursi.full + " : " + str(subj.visit) + " : CBF distributions" )
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path(os.path.join(perf_o_path, 'qc'))
            
            # move cbf native to subject's T1 space
            cbf_file = os.path.join(subj.derivatives_path('perf'), 'cbf', 'native_space', 'perfusion_calib.nii.gz')

            o_file_cbf = "VERIFY.cbf.T1.nii.gz"
            if not os.path.isfile(cbf_file):
                raise utils.CmdProcError("Missing: "+cbf_file)
                
            if not os.path.isfile(o_file_cbf):
                if space == 'TLRC':
                    
                    utils.Cmd("3dAllineate " +
                            " -master " + os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('T1w_SKSP.nii.gz')) +
                            " -1Dmatrix_apply " + os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('perf-to-T1.aff12.1D')) +
                            " -input " + cbf_file +
                            " -floatize" +
                            " -interp trilinear" +
                            " -final wsinc5"+
                            " -prefix " + o_file_cbf ).run()
                    
                else:
                    pCASLMbMd.apply_native_to_t1(subj, cbf_file, os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('T1w_brain.nii.gz')), o_file_cbf)
            
            
            
            # load wm mask
            i_file_wm = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-wm_mask.nii.gz'))
            if not os.path.isfile(i_file_wm):
                raise utils.CmdProcError("Missing: "+i_file_wm)
            img_wm = nib.load(i_file_wm)
            wm_mask = img_wm.get_data()
            
            # load gm mask
            i_file_gm = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-gm_mask.nii.gz'))
            if not os.path.isfile(i_file_wm):
                raise utils.CmdProcError("Missing: "+i_file_gm)
            img_gm = nib.load(i_file_gm)
            gm_mask = img_gm.get_data()

            # load T1 space CBF
            img_cbf = nib.load(o_file_cbf)
            cbf_data = img_cbf.get_data()

            # plot distributions
            
            # plot settings
            afniWarmCool = ["#f97306", "#056eee"]
            sns.set_palette(afniWarmCool)
        
            sns.set_style("white")
            sns.set_style("ticks")
            #plt.figure(figsize=(1.2, 1.0))
        
            plt.title("CBF masked distributions")
            plt.xlabel("CBF (ml/min/100ml)")
            # plot masked distributions
            sns.kdeplot( cbf_data[wm_mask.nonzero()].flatten(), shade=False, label="CBF in WM", linewidth=1.5 );
            sns.kdeplot( cbf_data[gm_mask.nonzero()].flatten(), shade=False, label="CBF in GM", linewidth=1.5 );
            # remove top/right axis lines
            sns.despine(top=True, right=True, left=True, bottom=False, offset=None, trim=False)
            plt.xlim(-50.0, 120.0)
            plt.legend()
            # display plots
            #plt.show()
            #plt.close()
            plt.savefig(subj.bids_file_join('CBF_in_mask.png'), dpi=None, facecolor='w', edgecolor='w', orientation='portrait', format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
            plt.close()
            
            utils.purge_paths(o_file_cbf)
        
        
        
        
        
        
        ### snapshot of spatial normalization
        def snapshot_spatial_normalization(subj, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
          
            try:
                # if not there, make it
                utils.make_path(os.path.join(subj.derivatives_path(bidsType='perf'), 'qc'))
                
                if not os.path.isfile(os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', subj.bids_file_join('perf', space+'.jpg'))):
                    
                    if space == 'MNI':
                        pCASLMbMd.apply_spatial_normalization(subj, 
                                                        os.path.join(subj.derivatives_path('perf'), 'topup.dewarped_fmap_Tmean.masked.nii.gz'),
                                                        os.path.join(subj.derivatives_path('perf'), 'qc', 'VERIFY.nii.gz'),
                                                        space='MNI', template_key='mni_icbm152_nlin_asym_2009c')
                        
                        ulay_file = utils.Env.template_lib[template_key]['head']
                        
                    else:
                        pCASLMbMd.apply_spatial_normalization(subj, 
                                                        os.path.join(subj.derivatives_path('perf'), 'topup.dewarped_fmap_Tmean.masked.nii.gz'), 
                                                        os.path.join(subj.derivatives_path('perf'), 'qc', 'VERIFY.nii.gz'), 
                                                        space='TLRC' )                
                        
                        ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/TT_N27+tlrc.HEAD"
                        
                        
                        
                    
                    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF : QC "+space+" Snapshots" )
                    olay_file = os.path.join(subj.derivatives_path('perf'), 'qc', 'VERIFY.nii.gz')
                    utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, olay_file, os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', subj.bids_file_join('perf', space)))).run()
                    
            except:
                raise
            finally:
                utils.purge_paths(os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', 'VERIFY*'))
            
        
        
        
        
        
        
        def aggregate_qc_pdf(subj, space='MNI'):
            
            """
            Collects QC images into a PDF with text
            """
            import pylatex
            from pylatex.utils import NoEscape        
            
            perf_path = subj.derivatives_path(bidsType='perf')
            
            o_file = os.path.join(perf_path, 'qc', 'qc-preprocessing')
            if not os.path.isfile(o_file + '.pdf'):
                
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF : QC PDF" )
                
                
                geometry_options = {"tmargin": "1cm", "lmargin": "3cm"}
                doc = pylatex.Document(page_numbers=False, geometry_options=geometry_options)
    
                doc.preamble.append(pylatex.Command('title', subj.study.label + ' PERF Preprocessing QC'))
                doc.preamble.append(pylatex.Command('author', subj.ursi.full + " : Visit " + str(subj.visit)))
                doc.preamble.append(pylatex.Command('date', NoEscape(r'\today')))
                doc.append(NoEscape(r'\maketitle'))
                doc.append(NoEscape(r'\newpage'))
                
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('Spatial Normalization', numbering=False)):
                        doc.append('Colored edge overlay should align with tissue boundaries')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image( os.path.join(perf_path, 'qc', subj.bids_file_join('perf_'+space+'.jpg')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('Distortion Corrected FMAP registration over '+space+' template')

                        
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('Distortion Correction', numbering=False)):
                        doc.append('Input images should contain minimal artifacts')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(perf_path, 'qc', 'topup_input_qc_onescl.sag.png'), width=NoEscape('0.6\\linewidth'), placement=NoEscape('\\centering') )
                            pic1.add_caption('Uncorrected FMAP Input - Row1-AP - Row2-PA - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(perf_path, 'qc', 'topup_output_qc_onescl.sag.png'), width=NoEscape('0.6\\linewidth'), placement=NoEscape('\\centering') )
                            pic2.add_caption('Corrected FMAP Output - Sagittal View')   
                
                
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('PERF Mask', numbering=False)):
                        doc.append('Mask should cover entire brain')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(perf_path, 'qc', 'mask.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('PERF Mask - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(perf_path, 'qc', 'mask.cor.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('PERF Mask - Coronal View')
                        with doc.create(pylatex.Figure(position='h')) as pic3:
                            pic3.add_image(os.path.join(perf_path, 'qc', 'mask.axi.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic3.add_caption('PERF Mask - Axial View')  
                        
                        
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('CBF', numbering=False)):
                        doc.append('Results')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(perf_path, 'qc', subj.bids_file_join('CBF.axi.png')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('CBF Native Space - Axial View')
                        doc.append(NoEscape(r'\newpage'))
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(perf_path, 'qc', subj.bids_file_join('CBF_in_mask.png')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('CBF by Tissue')   
                            
                            
                doc.generate_pdf(o_file, clean_tex=True, compiler='pdflatex')
    
    
    
    
    

class pCASL2d:

    def standard(subj, space='TLRC'):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ASL : Preprocessing" )
        
        check_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL-to-'+space+'.NWARP.nii.gz'))
        if not os.path.isfile(check_file):
            pCASL2d.motion_correction(subj)
            pCASL2d.extract_mask(subj)
            pCASL2d.extract_motion_estimates(subj)
            pCASL2d.spatial_normalization(subj, space=space)
        
        
        pCASL2d.generate_cbf(subj, min_tr_ms=4140, pld=1800)

        pCASL2d.qc.snapshot_mask(subj)
        pCASL2d.qc.cbf(subj)
        pCASL2d.qc.snapshot_cbf(subj)
        pCASL2d.qc.snapshot_spatial_normalization(subj, space=space)
        pCASL2d.qc.aggregate_qc_pdf(subj)

        pCASL2d.clean_proc(subj)
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ASL : Preprocessing Complete" )

            



    def motion_correction(subj):
        
        perf_raw_path = subj.bids_modality_path(bidsType='perf')
        if not os.path.isdir(perf_raw_path):
            raise ValueError("Missing input:\n"+perf_raw_path)
    
        perf_o_path = subj.derivatives_path(bidsType='perf')
        utils.make_path(perf_o_path)

        
        import shutil
        
        
        # get PD as M0
        src_file = os.path.join(perf_raw_path, subj.bids_file_join('acq-pcasl2d_m0scan.nii.gz'))
        M0_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.nii.gz'))
        shutil.copyfile(src_file, M0_file)
        
        # There is a bug in 2dImReg that causes it to misreport dimension errors with NIFTI input
        # Work around this by working in AFNI file format until we get beyond 2dImReg
        
        # Despike
        i_file = os.path.join(perf_raw_path, subj.bids_file_join('acq-pcasl2d_asl.nii.gz'))
        o_file = os.path.join(perf_o_path, subj.bids_file_join('dspk+orig.HEAD'))
        if not os.path.isfile(o_file):
            utils.prettyOut("Despike")
            utils.Cmd("3dDespike -NEW -ignore 0 -nomask -cut 4 4.5 -prefix " + o_file + " " + i_file).run()        

        
        # 2d slice registration
        i_file = o_file
        o_file = os.path.join(perf_o_path, subj.bids_file_join('2dreg.nii.gz'))
        if not os.path.isfile(o_file):
            utils.prettyOut("DATA IS REGISTERED IN 2D ")
            utils.Cmd("2dImReg  -input " + i_file + " -base 0 -prefix " + o_file + " -dprefix " + os.path.join(perf_o_path, subj.bids_file_join('pCASL.2dreg')) ).run()
            utils.purge_paths(i_file)
            
        
        # 3d registration to M0
        i_file = o_file
        o_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0-reg.nii.gz'))
        if not os.path.isfile(o_file):
            utils.prettyOut("DATA IS REGISTERED IN 3D ")
            utils.Cmd( "3dvolreg -twopass -twodup -verbose -prefix " + o_file + " -base " + M0_file + " " + i_file ).run()
            # convert back to short
            utils.Cmd("fslmaths {0} {0} -odt short".format(o_file)).run()


        
        


    def extract_mask(subj):
        #########################################################
        #     create brain mask and masked data
        #########################################################
        perf_o_path = subj.derivatives_path(bidsType='perf')
        mask_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.mask.nii.gz'))
        M0_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.nii.gz'))
        if not os.path.isfile(mask_file):
            utils.prettyOut("Base mask")
            utils.Cmd( "3dAutomask -q -clfrac 0.50 -peels 3 -nbhrs 15 -dilate 0 -prefix " + mask_file + " " + M0_file ).run()

    





    def extract_motion_estimates(subj):
        
        perf_raw_path = subj.bids_modality_path(bidsType='perf')
        perf_o_path = subj.derivatives_path(bidsType='perf')
        
        #########################################################
        #    estimate run motion
        #########################################################
        o_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.3dreg.log.raw'))
        p_raw_file = os.path.join(perf_raw_path, subj.bids_file_join('acq-pcasl2d_asl.nii.gz'))
        if not os.path.isfile(o_file):
            utils.prettyOut("Estimate run motion parameters ")
            utils.Cmd( "3dvolreg -dfile " + o_file + " -prefix NULL -base 0 " + p_raw_file ).run()


        #########################################################
        #    calculate FD
        #########################################################
        i_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.3dreg.log.raw'))
        o_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.3dreg.log.raw.deriv'))

        #    this can also be used as level1 motion regressor along with raw
        if not os.path.isfile(o_file):
            utils.prettyOut("calculate FD")
            utils.Cmd("1d_tool.py -infile "+i_file+" -derivative -write " + o_file).run()

        i_file = o_file
        o_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.3dreg.log.raw.deriv.fd.rowsum'))
        if not os.path.isfile(o_file):
            utils.Cmd("1deval -a "+i_file+"'[1]' -b "+i_file+"'[2]' -c "+i_file+"'[3]' -d "+i_file+"'[4]' -e "+i_file+"'[5]' -f "+i_file+"'[6]' -expr '100*sind(abs(a)/2) + 100*sind(abs(b)/2) + 100*sind(abs(c)/2) + abs(d) + abs(e) + abs(f)' > "+ o_file ).run()

        

        
        
    def align_to_T1(subj):
        """
        create alignment transform
        """
        
        import shutil
        
        try:
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            anat_path = subj.derivatives_path(bidsType='anat')
            tx_path = os.path.join(perf_o_path, subj.bids_file_join('pCASL-to-T1.aff12.1D'))
            
            if not os.path.isfile(tx_path):
    
                utils.prettyOut("create alignment transform")
                temp_path = os.path.join(perf_o_path, "regTemp_"+os.path.split(os.path.expanduser('~'))[-1])
                utils.make_path(temp_path)
    
                # save where we are now
                starting_path = os.getcwd()
                # move into that directoy
                os.chdir(temp_path)
    
                # copy visit T1 into that directory for the alignment
                shutil.copyfile(os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz')), 
                                os.path.join(temp_path, subj.bids_file_join('T1w_SKSP.nii.gz'))
                                )
                
                
                M0_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.nii.gz'))
                mask_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.mask.nii.gz'))
                pd_ma_file = os.path.join(temp_path, subj.bids_file_join('pCASL.M0.masked.nii.gz'))
                if not os.path.isfile(pd_ma_file):
                    utils.prettyOut("Base mask")
                    utils.Cmd( "fslmaths "+M0_file+" -mas "+ mask_file + " " + pd_ma_file).run()
                    
    
                utils.Cmd("align_epi_anat.py -epi2anat -anat "+os.path.join(temp_path, subj.bids_file_join('T1w_SKSP.nii.gz'))+" -epi "+ pd_ma_file + " -epi_base 0 -volreg off -epi_strip None -big_move -deoblique on -anat_has_skull no -tshift off -cost lpc -Allineate_opts -source_automask+4").run()
    
                # mv the tx file into perf
                shutil.copyfile( os.path.join(temp_path, subj.bids_file_join('pCASL.M0.masked_al_mat.aff12.1D')), tx_path )

                os.chdir(starting_path)
                utils.purge_paths(temp_path)

        except:
            raise
        




    def spatial_normalization(subj, space='MNI'):
        try:
            
            # make sure there is an affine transform
            pCASL2d.align_to_T1(subj)
            
            perf_o_path = subj.derivatives_path(bidsType='perf')
            anat_path = subj.derivatives_path(bidsType='anat')
            
            o_path = os.path.join(perf_o_path, subj.bids_file_join('pCASL-to-'+space+'.NWARP.nii.gz'))
            if not os.path.isfile(o_path):
                utils.Cmd("3dNwarpCat -warp1 {} -warp2 {} -warp3 {} -prefix {} ".format(os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.WARP.nii.gz')),
                                                                                        os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.aff12.1D')),
                                                                                        os.path.join(perf_o_path, subj.bids_file_join('pCASL-to-T1.aff12.1D')),
                                                                                        o_path)).run()
            
        except:
            raise        
        
        
        
        
        
    def apply_spatial_normalization(subj, i_path, o_path, master_path=None, interp='wsinc5', space='MNI'):
        """
        
        Parameters
        ----------
        subj : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        """
    
        try:
    
            if not os.path.isfile(o_path):
                
                perf_o_path = subj.derivatives_path(bidsType='perf')
                anat_path = subj.derivatives_path(bidsType='anat')
                
                if master_path == None:
                    master_path = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.'+space+'.nii.gz'))
                    if not os.path.isfile(master_path):
                        raise ValueError('apply_spatial_normalization: Could not find master file ' + master_path)

                utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(i_path,
                                                                                                       os.path.join(perf_o_path, subj.bids_file_join('pCASL-to-'+space+'.NWARP.nii.gz')),
                                                                                                       master_path,
                                                                                                       interp,
                                                                                                       o_path
                                                                                                       )).run()
            
        except:
            raise
    
    



        
    def blur(subj, fwhm=6):
        """
        blur in native space
        """

        i_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL.M0-reg.nii.gz'))
        o_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL.M0-reg.fwhm'+str(fwhm)+'.nii.gz'))
        if not os.path.isfile(o_file):
            utils.Cmd( "3dmerge -1blur_fwhm " + str(fwhm) + " -doall -prefix " + o_file + " " + i_file ).run()

        i_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL.M0.nii.gz'))
        o_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL.M0.fwhm'+str(fwhm)+'.nii.gz'))
        if not os.path.isfile(o_file):
            utils.Cmd( "3dmerge -1blur_fwhm " + str(fwhm) + " -doall -prefix " + o_file + " " + i_file ).run()
            
            
        
        
        
    def generate_cbf(subj, min_tr_ms=4140, pld=1800, T1_blood=1650):
    
        """
    
        20180730 - JL - The code for calculating CBF could be updating in the future to more faithfully follow the Alsop formula (2014)
        For the previous method see generate_pcasl_cbf_wang(). If you don't have a PD image you'll have to use generate_pcasl_cbf_wang()
        or change this function to loop over pairs so control image can be used as m0 normalizing image.
        
        
        Alsop et al. 2014:
    
              6000 * λ * (Control-Tagged) * e^(PLD/T1_blood)
        CBF = ---------------------------------------------------  [ml/100g/min]
              2 * alpha * T1_blood * PD * (1 - e^-(tau/T1_blood))
    
    
        Where:
            λ        (blood-brain partition coefficient)                 [0.9 mL/g]        
            PLD      (post-label delay)                                  [ms]
            T1_blood (longitudinal relaxation time of blood at 3.0 Tesla [1650 ms]
            α        (labeling efficiency) for PCASL                     [0.85]
            PD       (spin-free scaling image (M0))
            tau      (label duration of the labeling RF pulse train)     [# of RF blocks * block duration]
                       
            block duration = 0.0185s                                     [Function of the # RF pulses and gaps]
            
                
        
        The time between the end of this pulse train and image acquisition is referred to as the post labeling delay (PLD). 
        PLD refers to the time at which the end of the labeled bolus leaves the labeling plane in PCASL.
        For CBF quantification using PCASL, the ideal case is that the PLD is set just longer than
        the longest value of arterial transit time(ATT) present in the subject. Under these conditions, the entire labeled
        bolus is delivered to the tissue prior to image acquisition, and the CBF measurement will be
        unbiased by incomplete delivery. However, because the ASL signal decays with time
        constant T1 after labeling, it is too costly in terms of SNR to be extremely conservative in
        the choice of PLD such that PLD is guaranteed to be strictly longer than ATT under all
        circumstances.
        In healthy gray matter, ATT can vary between 500–1500ms depending on
        the labeling location and the tissue location in the brain, but in cerebrovascular disease and
        in deep white matter, ATT can be 2000ms or longer. The choice of PLD is therefore a
        compromise, such that SNR is acceptable, and that in the large majority of cases the ASL
        signal will accurately reflect CBF. However, it should be understood that areas of low ASL
        signal may reflect some combination of low CBF and unusually long ATT, and not specifically low CBF.
    
    
    
        brain/blood partition coefficient - This factor scales the signal intensity of tissue to that of blood.
        
        If TR is less than 5s, the PD image should be multiplied by the factor (1/(1 − e−TR/T1,tissue)), where T1,tissue is the
        assumed T1 of gray matter, in order to compensate for T1 relaxation. Using a reduced TR
        and T correction may potentially reduce errors associated with a brain-averaged λ.
    
         
    
    
    
        min_tr_ms
        Required
        See inline comments in code
    
        tag_delay_ms
        Required
        See inline comments in code
    
    
        """
    
        import os
        import numpy as np
        import numpy.ma as ma
        from termcolor import cprint    #    for colored printing
        import nibabel as nib
    
        try:
    
            ##################################################
            # OPEN IMAGES
            ##################################################
            perf_o_path = subj.derivatives_path(bidsType='perf')
            
            # fwhm = 6
            # pCASL2d.blur(subj, fwhm=6)            
            # asl_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0-reg.fwhm'+str(fwhm)+'.nii.gz'))
            # m0_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.fwhm'+str(fwhm)+'.nii.gz'))
            # cbf_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.CBF.fwhm'+str(fwhm)+'.nii.gz'))
            
            # no blur prior to CBF estimate
            asl_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0-reg.nii.gz'))
            m0_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.nii.gz'))
            cbf_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.CBF.nii.gz'))
            

            if not os.path.isfile(cbf_file):
                    
                print("Loading data...")
                mask_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.mask.nii.gz'))
                mask_img = nib.load(mask_file)
                mask_data = mask_img.get_fdata().astype('int16')
                
                # load pcasl volume
                asl_img = nib.load(asl_file)
                # unpack dimensions to make code more readable
                (asl_xdim,asl_ydim,asl_zdim,asl_tdim) = asl_img.get_fdata().shape
                asl_data = np.asarray(asl_img.dataobj).astype('int16')
                
                # load perf volume
                m0_img = nib.load(m0_file)
                m0_data = m0_img.get_fdata().astype('int16')
                # when i create m0 vol i'm do so as a masked array so that we can use it in division.
                m0_vol = ma.masked_where(mask_data == 0, m0_data, copy=True)

                
                
                ##################################################
                # DEFINE VALUES FOR CBF CALCULATION
                ##################################################
        
                # blood/tissue water partition coefficient
                # 0.9 for pCASL Alsop 2014
                # [mL/g]
                blood_tissue_h2o_part_coef = 0.9
        
                # keep clear what these are in calculations
                # Number of RF Blocks = "Num RF Blocks" from image time dimension
                num_rf_blocks = asl_tdim
                
                # number of slices from image z dimension
                num_slices = asl_zdim
        
                # rf_block_duration
                # [seconds]
                rf_block_duration = 0.0185
        
                # label_duration (aka Labeltime) =  Number of RF Blocks*rf_block_duration;
                label_duration = (num_rf_blocks*(rf_block_duration*1000))
                
                # Post tagging delay/Post-labeling delay; get from sequence protocol
                # [msec]
                post_tag_delay_time = pld
        
                # Slicetime
                # time for getting one slice, which can calculated by
                # phase encoding lines*echo spacing time [+ fat saturation time and crusher time+ slice selection gradient time + phase refocusing gradient time]
                # RORDEN: "To compute slicetime get the minimal TR by clicking the "TR" window in the protocol panel, then slicetime=(minTR-labelingtime-delaytime)/#slices.
                # For example, the MCBI default sequence has a minimum TR= 2090ms plus the delay time (e.g. with a 1200ms delay, the minimum time is 3290ms)
                # For example if MinTR=3090, labelingtime=1480, delaytime=1000, slices=17, then slicetime= 35.88235294"
        
                # minTR               = get the minimal TR by clicking the "TR" window in the protocol panel
                # tagging_time        = see above
                # post_tag_delay_time = see above
                #
                # Slicetime = Slice acquisition time:msec
                # [msec]
                slice_time = (min_tr_ms-label_duration-post_tag_delay_time)/float(num_slices)
        
                # r1a=1/BloodT1; %0.601sec for 3T
                # the longitudinal relaxation time of blood
                # 1650 for pCASL Alsop 2014
                # [msec]
                #T1_blood = 1650
                
                # labeling efficiency, 0.85 for pCASL Alsop 2014
                labeling_efficiency = 0.85
        
        
                print("\n")
                cprint("\n=======================================================", 'red', 'on_grey', attrs=['bold'] )
                cprint("  CBF Parameters", 'red', 'on_grey', attrs=['bold'] )
                cprint('  + blood_tissue_h2o_part_coef : {0:.2f}'.format(blood_tissue_h2o_part_coef), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + rf_block_duration : {0:.6f}'.format(rf_block_duration), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + pld : {0:d}'.format(pld), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + min_tr_ms : {0:d}'.format(min_tr_ms), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + num_rf_blocks : {0:d}'.format(num_rf_blocks), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + num_slices : {0:d}'.format(num_slices), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + label_duration : {0:.2f}'.format(label_duration), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + slice_time : {0:.2f}'.format(slice_time), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + T1_blood : {0:.2f}'.format(T1_blood), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + labeling_efficiency : {0:.2f}'.format(labeling_efficiency), 'red', 'on_grey', attrs=['bold'] )
                cprint('  + using PD image for m0 normalization', 'red', 'on_grey', attrs=['bold'] )
                cprint("=======================================================", 'red', 'on_grey', attrs=['bold'] )

        
                #    we are creating a slice timing adjustment for the total delay
                #    which includes the post-tagging delay plus the timing of each slice acquistion from bottom to top of the brain
                #    don't use with slice timing corrected input volumes [EITHER/OR]. This method used for ASLToolbox comparison
                #
                #    create a slicetimearray of ones (64*64,24)
                #    also add delay time aka post_label delay
                slice_time_array = np.ones(asl_xdim*asl_ydim*asl_zdim,dtype=float).reshape(asl_xdim,asl_ydim,asl_zdim)
                slice_weight = (np.arange(0,asl_zdim,1)*slice_time + post_tag_delay_time)
                for slice_indx in range(asl_zdim):
                    slice_time_array[:,:,slice_indx] = slice_weight[slice_indx]
                
        
                # perfusion image is mean of all labeled images subtracted from all control or unlabeled images (control - label)
                mean_diff = np.mean(asl_data[:,:,:,1::2], axis=3) - np.mean(asl_data[:,:,:,::2], axis=3)
        
                # calculate mean cbf using alsop formula
                #          6000 * λ * (Control-Tagged) * e^(PLD/T1_blood)
                #    CBF = ---------------------------------------------------  [ml/100g/min]
                #          2 * alpha * T1_blood * PD * (1 - e^-(tau/T1_blood))        
        
                # The factor of 6000 converts the units from ml/g/s to ml/(100g)/min, which is customary in the physiological literature.
                # The 1000 converts from milleseconds to seconds. The input parameters are given in milleseconds in Alsop 2014.
                # To help in reporting i've left the inputs in milleseconds and the two rescaling factors separate in the formula. i.e. inputs in milleseconds, output in minutes
        
                # calculate mean CBF: no slice time correction
                #mean_cbf_vol = ( 6000 * 1000 * blood_tissue_h2o_part_coef * mean_diff * np.exp(pld/T1_blood) ) / ( 2 * labeling_efficiency * T1_blood * m0_vol * (1-(np.exp(-(label_duration/T1_blood)))) )
                
                # divide by zero errors could be a good diagnostic but i have found it to be driven by non-parenchyma while steps below clean up the debris
                # therefore, turn off these errors for the calculation
                with np.errstate(divide='ignore', invalid='ignore'):
                    # calculate mean CBF: includes slice time correction
                    mean_cbf_vol = ( 6000 * 1000 * blood_tissue_h2o_part_coef * mean_diff * np.exp( (slice_time_array)/T1_blood) ) / ( 2 * labeling_efficiency * T1_blood * m0_vol * (1-(np.exp(-(label_duration/T1_blood)))) )
                
                # convert INF to 0 to forestall conversion to +/-1.7976931348623157e+308 by np.nan_to_num().
                mean_cbf_vol[np.isinf(mean_cbf_vol)] = 0.0
                
                # remove nan's that might occur above from divide-by-zero
                # this will also convert INF to +/-1.7976931348623157e+308 which i don't find useful as it won't store in float32
                # and is a useless arfifact of divide by zero
                mean_cbf_vol = np.nan_to_num(mean_cbf_vol).astype('float32')

                # finally, write out mean cbf volume
                print("Writing Mean CBF Volume...")
                masked_data = mean_cbf_vol.data * mask_data
                # remove negative CBF which makes no sense and implausibly high CBF. oxford_asl also zeros out negative values
                masked_filtered_data = np.where( (masked_data < 0.0) | (masked_data > 200.0), 0, masked_data)
                # create nifti image object
                o_img = nib.nifti1.Nifti1Image(masked_filtered_data, m0_img.affine, m0_img.header)
                # save the image to file
                nib.save(o_img, cbf_file)

                
        except:
    
            raise
    
    


    def clean_proc(subj):
        """
        
        Parameters
        ----------
        subj : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        """
    
        try:
            
            utils.purge_paths( os.path.join(subj.derivatives_path('perf'), subj.bids_file_join('dspk+orig*')),
                               os.path.join(subj.derivatives_path('perf'), subj.bids_file_join('2dreg*'))
                              )
                                  
            
        except:
            raise





    class qc:
        
        
        ### snapshot of spatial normalization
        def snapshot_spatial_normalization(subj, space='MNI'):
            """
            """    
            
            try:
                if not os.path.isfile(os.path.join(subj.derivatives_path('perf'), 'qc', subj.bids_file_join('perf_'+space+'.jpg'))):
                    
                    utils.make_path( os.path.join(subj.derivatives_path('perf'), 'qc') )
                    
                    verify_file = os.path.join(subj.derivatives_path('perf'), 'qc', 'VERIFY.nii.gz')
                    if not os.path.isfile(verify_file):
                        pCASL2d.apply_spatial_normalization(subj, 
                                                        os.path.join(subj.derivatives_path('perf'), subj.bids_file_join('pCASL.M0.nii.gz')), 
                                                        verify_file, space=space)
                
            
                        
                    if space == 'MNI':
                        ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/AFNI_MNI152_T1_2009c+tlrc.HEAD"
                    else:
                        ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/TT_N27+tlrc.HEAD"
                        
                        
                    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF : QC "+space+" Snapshots" )
        
                    utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, verify_file, os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', subj.bids_file_join('perf', space)))).run()
                    
            except:
                raise
                
            finally:
                utils.purge_paths(os.path.join(subj.derivatives_path(bidsType='perf'), 'qc', 'VERIFY*'))
                
            
            
            
            
        def snapshot_mask(subj):
            """
            create mask images into qc directory
            """
            perf_o_path = subj.derivatives_path(bidsType='perf')
            utils.make_path( os.path.join(perf_o_path, 'qc') )
            
            if not os.path.isfile(os.path.join(perf_o_path, 'qc', subj.bids_file_join('pCASL_mask.axi.png'))):
                
                ulay_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.M0.nii.gz'))
                olay_file = os.path.join(perf_o_path, subj.bids_file_join('pCASL.mask.nii.gz'))
    
                utils.make_path(os.path.join(perf_o_path, 'qc'))
                
                utils.Cmd("@chauffeur_afni \
                        -ulay "+ulay_file+" \
                        -olay "+olay_file+" \
                        -prefix "+ os.path.join(perf_o_path, 'qc', 'mask') +" \
                        -set_xhairs OFF \
                        -montx 4 \
                        -monty 1 \
                        -delta_slices 4 4 2").run()

            
            
            
        def snapshot_cbf(subj):
            
            
            """
            
            """    
            
            if not os.path.isfile(os.path.join(subj.derivatives_path('perf'), 'qc', subj.bids_file_join('CBF.axi.png'))):
                
                ulay_file = os.path.join(subj.derivatives_path('perf'), subj.bids_file_join('pCASL.M0.nii.gz'))
                olay_file = os.path.join(subj.derivatives_path('perf'), subj.bids_file_join('pCASL.CBF.nii.gz'))
                utils.make_path(os.path.join(subj.derivatives_path('perf'), 'qc'))
                
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ASL : QC CBF Snapshots" )
                
                utils.Cmd("@chauffeur_afni \
                        -ulay "+ulay_file+" \
                        -olay "+olay_file+" \
                        -prefix "+ os.path.join(subj.derivatives_path('perf'), 'qc', subj.bids_file_join('CBF')) +" \
                        -set_xhairs OFF \
                        -montx 6 \
                        -monty 3 \
                        -delta_slices 1 1 1 \
                        -opacity 9 \
                        -cbar_topval 80 \
                        -cbar_ncolors 9 \
                        -func_range 60 -thr_olay 0.1").run()
                    


            
            
        def cbf(subj, fwhm=6):
            
            """
            ALSOP 2015
            https://doi.org/10.1002/mrm.25197
        
            For evaluating the quality of ASL MRI images in clinical practice we advise the following
            checks:
            
            1.
            For PCASL scans, look for areas of low labeling efficiency. First, identify which
            arteries should have been labeled. Typically, this will include internal and external
            carotid arteries, and vertebral arteries. If an angiogram is available, this can be used
            to verify the list of labeled arteries. Checking the Circle of Willis anatomy may also
            be of use in matching vascular territories to labeled arteries. When the labeling
            efficiency is low in an artery, the entire flow territory of that artery will
            demonstrate a low calculated CBF. When a low CBF area is seen that matches an
            entire vascular territory, with no apparent compensation from other arteries, a
            labeling failure should be considered, though this does not preclude the possibility
            of truly low CBF conditions, or abnormally long ATT. Labeling failures can be
            caused by tortuous vessels or resonance offsets in the labeling plane. The former
            may be addressed by adjusting the location of the labeling plane, in which case an
            additional angiogram would be helpful. The latter is commonly caused by dental
            work, and may be suggested by signal dropouts around the teeth in other images
            from this patient. Methods to address resonance offset related labeling problems in
            PCASL are discussed above.
            
            2.
            Note the overall gray matter CBF value. Absolute CBF values obtained in gray
            matter can vary significantly, even among healthy young adults, due to natural
            inter-subject and intra-subject variations. In addition, average numbers are sensitive
            to partial volume effects and the methods used for isolating gray matter signal. As a
            general rule, gray matter CBF values from 40–100 ml/min/100ml can be normal.
            When the overall gray matter CBF value is inconsistent with the expected values
            for the patient population, consider the possibility that there is a global reduction in
            labeling efficiency, or that the PD scan used for normalization was incorrectly
            acquired or scaled. Clear contrast between gray and white matter should be present,
            and if not, may signify either poor labeling or motion artifacts.
            
            3.
            Check for motion artifacts. As a subtractive technique, ASL is motion sensitive,
            though this sensitivity is mitigated by BS as discussed above. The presence of
            signal outside of the brain, frequently recognizable as signal from layers of skin or
            fat is a clear indication of significant subject motion. When possible, it may be
            useful to check individual label/control difference images before averaging to see
            whether artifacts arise from only a minority of these difference images. If so, these
            images can be excluded from the CBF calculation. In addition, motion correction
            by means of automated image registration algorithms can be performed, though
            these may not be effective when BS is very efficient, or when applied to the label/
            control difference images, as the individual image SNR in these cases is low. When
            BS is not used or is incomplete, image registration is likely to be more effective,
            but BS is nevertheless recommended as a primary means of reducing physiological
            noise and motion artifacts. In the ideal case, prospective motion correction methods
            can be used when available to reduce motion artifacts during acquisition (73), and
            some of these methods are compatible with background suppression.
            
            4.
            Look for intravascular artifacts. Hyperintense spots and serpiginous regions often
            represent intravascular signal. When observed, it is advisable to verify that the PLD
            was appropriate for the patient (see Table 1), as a low PLD will naturally generate
            ASL signals in larger arteries. Intra-arterial signal with a correct PLD suggests that
            delivery of labeled blood to tissue is delayed, through slow flow and/or circuitous
            or collateral routes of circulation. Intra-venous ASL signal suggests that an
            arteriovenous shunt is present. Note that CBF calculations over whole brain or
            large regions of interest may still be valid in the presence of intravascular artifact as
            long as flow crushing gradients were not used.
            
            5. 
            Check the borderzone (watershed) regions. The borderzone or watershed areas are
            at the more distal portions of each vascular territory, and will naturally have a
            longer ATT than other portions of the territory. Note that it is possible for low ASL
            signal in these regions to represent long ATT rather than low CBF, and an
            additional scan with longer PLD may help to distinguish between these two
            possibilities. An example of this effect is shown in Figure 9.
            """
        
            import nibabel as nib
            import os
            import seaborn as sns
            import matplotlib.pyplot as plt
        
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : QC CBF distributions" )
        
            if not os.path.isfile(os.path.join(subj.derivatives_path('perf'), 'qc', subj.bids_file_join('CBF_in_mask.png'))):
                    
                utils.make_path(os.path.join(subj.derivatives_path('perf'), 'qc'))
                os.chdir( os.path.join(subj.derivatives_path(bidsType='perf'), 'qc') )
                
                # move cbf native to subject's T1 space
                cbf_file = os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL.CBF.nii.gz'))
                o_file_cbf = "VERIFY.cbf.T1.nii.gz"
                if not os.path.isfile(cbf_file):
                    raise utils.CmdProcError("Missing: "+cbf_file)
                    
                if not os.path.isfile(o_file_cbf):
                    utils.Cmd("3dAllineate " +
                            " -master " + os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('T1w_SKSP.nii.gz')) +
                            " -1Dmatrix_apply " + os.path.join(subj.derivatives_path(bidsType='perf'), subj.bids_file_join('pCASL-to-T1.aff12.1D')) +
                            " -input " + cbf_file +
                            " -floatize" +
                            " -interp trilinear" +
                            " -final wsinc5"+
                            " -prefix " + o_file_cbf ).run()
            
                
    
                # load wm mask
                i_file_wm = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-wm_mask.nii.gz'))
                if not os.path.isfile(i_file_wm):
                    raise utils.CmdProcError("Missing: "+i_file_wm)
                img_wm = nib.load(i_file_wm)
                wm_mask = img_wm.get_data()
                
                # load gm mask
                i_file_gm = os.path.join(subj.derivatives_path(bidsType='anat'), subj.bids_file_join('spm-gm_mask.nii.gz'))
                if not os.path.isfile(i_file_wm):
                    raise utils.CmdProcError("Missing: "+i_file_gm)
                img_gm = nib.load(i_file_gm)
                gm_mask = img_gm.get_data()
                
                # load T1 space CBF
                img_cbf = nib.load(o_file_cbf)
                cbf_data = img_cbf.get_data()
                
                # plot distributions
                
                # plot settings
                afniWarmCool = ["#f97306", "#056eee"]
                sns.set_palette(afniWarmCool)
            
                sns.set_style("white")
                sns.set_style("ticks")
                #plt.figure(figsize=(1.2, 1.0))
            
                plt.title("CBF masked distributions")
                plt.xlabel("CBF (ml/min/100ml)")
                # plot masked distributions
                sns.kdeplot( cbf_data[wm_mask.nonzero()].flatten(), shade=False, label="CBF in WM", linewidth=1.5 );
                sns.kdeplot( cbf_data[gm_mask.nonzero()].flatten(), shade=False, label="CBF in GM", linewidth=1.5 );
                # remove top/right axis lines
                sns.despine(top=True, right=True, left=True, bottom=False, offset=None, trim=False)
                plt.xlim(-50.0, 120.0)
                plt.legend()
                # display plots
                #plt.show()
                #plt.close()
                plt.savefig(subj.bids_file_join('CBF_in_mask.png'), dpi=None, facecolor='w', edgecolor='w', orientation='portrait', format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
                plt.close()
                utils.purge_paths(o_file_cbf)
                
            
            
            
            
            
        def aggregate_qc_pdf(subj, space='TLRC'):
            
            """
            Collects QC images into a PDF with text
            """
            import pylatex
            from pylatex.utils import NoEscape        
            
            perf_path = subj.derivatives_path(bidsType='perf')
            
            o_file = os.path.join(perf_path, 'qc', 'qc-preprocessing')
            if not os.path.isfile(o_file + '.pdf'):
                
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : PERF : QC PDF" )
                
                
                geometry_options = {"tmargin": "1cm", "lmargin": "3cm"}
                doc = pylatex.Document(page_numbers=False, geometry_options=geometry_options)
    
                doc.preamble.append(pylatex.Command('title', subj.study.label + ' PERF Preprocessing QC'))
                doc.preamble.append(pylatex.Command('author', subj.ursi.full + " : Visit " + str(subj.visit)))
                doc.preamble.append(pylatex.Command('date', NoEscape(r'\today')))
                doc.append(NoEscape(r'\maketitle'))
                doc.append(NoEscape(r'\newpage'))
                
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('Spatial Normalization', numbering=False)):
                        doc.append('Colored edge overlay should align with tissue boundaries')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image( os.path.join(perf_path, 'qc', subj.bids_file_join('perf_'+space+'.jpg')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('Distortion Corrected FMAP registration over '+space+' template')
                
                
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('PERF Mask', numbering=False)):
                        doc.append('Mask should cover entire brain')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(perf_path, 'qc', 'mask.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('PERF Mask - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(perf_path, 'qc', 'mask.cor.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('PERF Mask - Coronal View')
                        with doc.create(pylatex.Figure(position='h')) as pic3:
                            pic3.add_image(os.path.join(perf_path, 'qc', 'mask.axi.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic3.add_caption('PERF Mask - Axial View')  
                        
                        
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('CBF', numbering=False)):
                        doc.append('Results')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(perf_path, 'qc', subj.bids_file_join('CBF.axi.png')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('CBF Native Space - Axial View')
                        doc.append(NoEscape(r'\newpage'))
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(perf_path, 'qc', subj.bids_file_join('CBF_in_mask.png')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('CBF by Tissue')   
                            
                            
                doc.generate_pdf(o_file, clean_tex=True, compiler='pdflatex')
  
        
  


        


        
        
        def plots(subj, mask_path):
            
            """
            Assumes data and mask in tlrc space so no transforms
            """
        
            from mayerlab import utils
            import os
            import seaborn as sns
            import matplotlib.pyplot as plt
        
            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : CBF distributions" )
        
        
            # load mask
            if not os.path.isfile(mask_path):
                raise utils.CmdProcError("Missing: "+mask_path)
            img_mask = utils.AfniImage.load(mask_path)
            mask = img_mask.data[:,:,0:38]
        
            
            # load T1 space CBF
            i_file_cbf = os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.mean.cbf+tlrc.HEAD")
            img_cbf = utils.AfniImage.load(i_file_cbf)
            cbf_data = img_cbf.data[:,:,0:38]
            utils.Cmd('3dBrickStat -slow -nan -count -mask ' + mask_path + ' ' + i_file_cbf).run()
            utils.Cmd('3dBrickStat -min -max -mask ' + mask_path + ' ' + i_file_cbf).run()
            
            # plot distributions
            
            # plot settings
            afniWarmCool = ["#f97306", "#056eee"]
            sns.set_palette(afniWarmCool)
        
            sns.set_style("white")
            sns.set_style("ticks")
            #plt.figure(figsize=(1.2, 1.0))
        
            plt.title("CBF masked distributions: " + subj.ursi.short)
            plt.xlabel("CBF (ml/min/100ml)")
            # plot masked distributions
            sns.kdeplot( cbf_data[mask.nonzero()].flatten(), shade=False, label="CBF in mask", linewidth=1.5 );
            # remove top/right axis lines
            sns.despine(top=True, right=True, left=True, bottom=False, offset=None, trim=False)
            #plt.xlim(-50.0, 120.0)
            # display plots
            plt.show()

    
    
    
    
        
        
    
        

def nd_blur(d, fwhm):

    """
    Gaussian 3d blur of numpy array
    JL: This technique gives results that exactly match what comes out of 3dmerge -1blur_fwhm
    for anat volume

    WIP

    """

    import numpy as np
    from scipy import ndimage
    import math

    if not isinstance(d,np.ndarray):
        raise ValueError('"d" not numpy array')

    if not isinstance(fwhm, (int, float)):
        raise ValueError('fwhm not correct format')

    # calculate sigma from fwhm. Filter below expects sigma for kernel width
    sigma = fwhm/math.sqrt(8*math.log(2))

    #delta = img.header.get_zooms()
    delta = (3.4375, 3.4375002, 6.0)
    delta = (3.43, 3.43, 6.00)
    sigma_seq = [(fwhm/z)/math.sqrt(8*math.log(2)) for z in delta]
    df_gb = ndimage.gaussian_filter(d,sigma=[0.6176885820276502, 0.6176885820276502, 0.35388408345334127,0.00])
    df_gb = ndimage.gaussian_filter(d,sigma=[0.6176885820276502, 0.6176885820276502])

    # if array dtype is int convert to float for precision
    if d.dtype.kind == 'i':

        # convert image dtype before blurring
        df = d.astype(np.float32)
        # apply filter
        df_gb = ndimage.gaussian_filter(df, sigma=sigma)
        # round away the decimals
        df_gb_r = df_gb.round(decimals=0)
        # convert back to int16
        d_write = df_gb_r.astype(d.dtype)

    # if array is float, proceed without dtype conversion
    elif d.dtype.kind == 'f':

        # blur
        d_write = ndimage.gaussian_filter(d, sigma=sigma)

    else:

        raise ValueError('d contains unknown data type')

    return d_write

