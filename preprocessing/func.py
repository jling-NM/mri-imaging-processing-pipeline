#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@package func
Preprocessing module for diffusion imaging
Author: jling@mrn.org

Created on Mon Mar 29 16:28:37 2021

@author: jling
"""


from mayerlab.preprocessing import utils
import os




class pipelines:
    
    def standard(subj):
        
        """
        
    
        Parameters
        ----------
        subj : TYPE
            fully process each task.
            
            despike(subj,study,run)
            timeshift(subj,study,run)
            imReg_2d(subj,study,run)
            imReg_3d(subj,study,run)
            estimate_run_motion(subj,study,run)
                
            #    finalize each run type
            #    condense each run into a final preproc file
            #    ready for level1
            #
            #    combine runs into single volume
            #    apply field distortion correction to unwarp data
            #    mask for run_type (from SBRef)
            #    affine alignment. Keep separate so that user can redo just this step
            #    full concatenated alignment path through affine and non-linear
            
            
        Returns
        -------
        None.
    
        """
            
        try:
            
            # all get tasks and runs from protocol
            proc_map = subj.get_subj_task_map()
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing" )
            
            for task, run_number in proc_map.items():

                if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
                    
                    if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc.nii.gz'))):
    
                        # do this before each task is run
                        # pull sbref for this task
                        pull_sbref(subj, task)
                
                        for run in range(1, run_number + 1):
                            # do this for each run            
                            despike(subj, task, run)
                            timeshift(subj, task, run)
                            imReg_2d(subj, task, run)
                            imReg_3d(subj, task, run)
                            estimate_run_motion(subj, task, run)
                   
                    
                    # do this for the task after each run processed
                    concate_runs(subj, task, run_list=range(1, run_number + 1))
                    create_motion_regressors(subj, task, run_list=range(1, run_number + 1))
                    create_motion_deriv_file(subj, task, run_list=range(1, run_number + 1))
                    estimate_field_distortion(subj, task)
                    apply_distortion_correction_sbref(subj, task)
                    generate_task_mask(subj, task)
                    apply_distortion_correction_concat(subj, task)
                    qc.snapshot_mask(subj, task)
                    clean_proc(subj, task)
                    

                # for each task regardless of 'preproc-dc.nii.gz' status
                spatial_normalization(subj, task)
                qc.snapshot_spatial_normalization(subj, task)
                    
                
                
            # after all tasks complete
            qc.aggregate_qc_pdf(subj)
            qc.run_mriqc(subj)
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing Complete" )
                
            
        except:
            raise
            
            
            
            
            
        
    def laptop(subj):
        try:
            
            # for LAPTOP you need the whole process
            # get tasks and runs from protocol
            proc_map = subj.get_subj_task_map()
    
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing" )
            
            
            for task, run_number in proc_map.items():
                
                if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
                    
                    if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc.nii.gz'))):
                        
                        # do this before each task is run
                        # pull sbref for this task
                        pull_sbref(subj, task)
                
                        for run in range(1, run_number + 1):
                            # do this for each run            
                            despike(subj, task, run)
                            timeshift(subj, task, run)
                            imReg_2d(subj, task, run)
                            imReg_3d(subj, task, run)
                            estimate_run_motion(subj, task, run)
                   
                    
                    # do this for the task after each run processed
                    concate_runs(subj, task, run_list=range(1, run_number + 1))
                    create_motion_regressors(subj, task, run_list=range(1, run_number + 1))
                    create_motion_deriv_file(subj, task, run_list=range(1, run_number + 1))
                    estimate_field_distortion(subj, task)
                    apply_distortion_correction_sbref(subj, task)
                    generate_task_mask(subj, task)
                    apply_distortion_correction_concat(subj, task)
                    spatial_normalization(subj, task, space='TLRC')
                    qc.snapshot_mask(subj, task)
                    qc.snapshot_spatial_normalization(subj, task, space='TLRC')
                    qc.run_mriqc(subj)
                    
                    clean_proc(subj, task)                    
                
                    
            qc.aggregate_qc_pdf(subj)
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing Complete" )
        except:
            raise        






def pull_sbref(subj, task):
    """
    

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    import shutil
    func_path = subj.derivatives_path('func')
    utils.make_path(func_path)
    
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='sbref', task=task, extension='nii.gz')
    if not os.path.isfile(i_file):
        raise ValueError("Cannot pull" + i_file)
        
    o_file = os.path.join(func_path, os.path.basename(i_file))
    if not os.path.isfile(o_file):
        shutil.copyfile(i_file, o_file)

    
    
    
    
def slice_timing_to_file(subj, i_file_json, o_file_txt):
    """
    Extract slice timeing from bids sidecar into text file

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    import json
        
    try:
        with open(i_file_json, 'r') as json_file:
            side_car_data = json.load(json_file)
        
        with open(o_file_txt, 'w') as timing_file:
            timing_file.write(" ".join(map(str, side_car_data['SliceTiming'])))
    
        return 
    except:
        raise
        
        
        
        
def despike(subj, task, run, ignore=1):
    """
    despike

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    ignore : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')
    utils.make_path(func_path)
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='nii.gz')
    o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
    
    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Despike")
        utils.Cmd("3dDespike -NEW -ignore " + str(ignore) + " -nomask -cut 4 4.5 -prefix "+o_file+" " + i_file).run()
        
        
        

def timeshift(subj, task, run, ignore=1):
    """
    timeshift

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    ignore : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')  
    
    i_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
    o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tshft.nii.gz'))

    timing_json = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='json')
    timing_txt = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'slice_timing.1D'))
    slice_timing_to_file(subj, timing_json, timing_txt)

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Slice Timing Correction")
        utils.Cmd("3dTshift -tpattern @"+timing_txt+" -ignore "+ str(ignore) + " -tzero 0 -prefix "+ o_file +" " + i_file).run()




def imReg_2d(subj, task, run):
    """
    2d registration to SBREF base file

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    reg_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg.nii.gz'))
    reg_params_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg'))

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : 2D Registration")

        # weird workaround for tshft output to avoid dataset dimension error
        utils.Cmd("3dcopy -overwrite " + base_file + " " + os.path.join(subj.derivatives_path('func'), 'tmp_'+task+'_sbref')).run()
        utils.Cmd("3dcopy -overwrite " + reg_file + " " + os.path.join(subj.derivatives_path('func'), 'tmp_'+task+'_dspk')).run()
        #utils.Cmd("2dImReg -input " + reg_file + " -basefile " + base_file + " -base 0 -prefix " + o_file + " -dprefix " + reg_params_file).run()
        utils.Cmd("2dImReg -input " + os.path.join(subj.derivatives_path('func'), 'tmp_'+task+'_dspk+orig') + " -basefile " + os.path.join(subj.derivatives_path('func'), 'tmp_'+task+'_sbref+orig') + " -base 0 -prefix " + o_file + " -dprefix " + reg_params_file).run()
        utils.purge_paths(os.path.join(subj.derivatives_path('func'), 'tmp_'+task+'*'))
        


        
        
def imReg_3d(subj, task, run):
    """
    3d registration of processed data to SBREF base file to correct motion

    This outputs parameters used as motion regressors. See create_run_motion_regressors()
    ".Movement.Regressor.mc.1D"):

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreg.nii.gz'))
    mat_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'mat.mc'))
    mot_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'log.mc'))
    
    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : 3D Registration")
        utils.Cmd("3dvolreg -dfile " + mot_file + " -1Dmatrix_save " + mat_file + " -prefix " + o_file + " -base " + base_file + " " + i_file ).run()




def estimate_run_motion(subj, task, run):
    """
    Estimate movement parameters
    Register the raw, unprocessed run data to SBRef to estimate motion.

    This is used for motion quantification. See create_run_motion_regressors()
    ".Movement.Regressor.1D"):

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='nii.gz')
    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreglog.raw'))

    #    only run if output does not already exist
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + ' run-' + str(run) + " : Estimate motion parameters")
        utils.Cmd("3dvolreg -dfile " + o_file + " -prefix NULL -base 0 -base " + base_file + " " + i_file).run()





def concate_runs(subj, task, run_list, ignore=1):
    """
    

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')
    run_input_paths = []
    for run_number in run_list:
        run_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task, 'run-'+str(run_number), '3dreg.nii.gz')))

    o_file = os.path.join(func_path, subj.bids_file_join('task-' + task, 'preproc-uncorrected.nii.gz'))

    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : Concatenate Runs")
        utils.Cmd("3dTcat -prefix " + o_file + " " + " ".join([file + "'["+str(ignore)+"-$]'" for file in run_input_paths])).run()
        
        
        
    
    
    
def create_motion_regressors(subj, task, run_list, ignore=1):
    """
    creates file appropriate for level1

    """

    mot_quant_input_paths = []
    mot_reg_input_paths = []
    for run_number in run_list:
        mot_quant_input_paths.append(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), '3dreglog.raw')))
        mot_reg_input_paths.append(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), 'log.mc')))
        
    # original movement regressor created with raw epi input. This is used for motion quantification
    # create from available run data
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor.1D'))
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor.1D")
        utils.Cmd("tail -q -n +" + str(ignore+1) + " " + ' '.join(mot_quant_input_paths) +" > " + o_file).run()


    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc.1D'))
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor-mc.1D")
        utils.Cmd("tail -q -n +" + str(ignore+1) + " " + ' '.join(mot_reg_input_paths) +" >> " + o_file).run()


    

    
def create_motion_deriv_file(subj, task, run_list):

    """
    Create derivative of motion regressor file

    """

    i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc.1D'))    
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc-deriv.1D'))

    run_length = []
    for run_number in run_list:
        with open(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), '3dreglog.raw'))) as f:
            run_length.append(str(sum(1 for _ in f) - 1))
            

        
    #    generate motion derivative regessors if they don't exist
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor-mc-deriv.1D")
        utils.Cmd("1d_tool.py"
                " -infile " + i_file +
                " -set_run_lengths " +  ' '.join(run_length) +
                " -derivative"
                " -write " + o_file).run()
        
        
        
        
        
        
def estimate_field_distortion(subj, task):
    """
    Estimate distortion from fmaps

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # locate fmaps for task. if there is no fmap set specifically for task, use 'fmri'
    # approach is to give each task its own map set even if it is just duplication.
    # perhaps a better approach would be to map these in the visit settings or the BIDS sidecar attributed intended for such
    #
    
    
    # if there is neither a complete set of task specific or fmri general fmaps, issue a warning then skip distoration estimate
    ap_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-'+task+'_dir-ap_epi.nii.gz') )
    if not os.path.isfile(ap_file):
        ap_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-fmri_dir-ap_epi.nii.gz') )
        
    pa_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-'+task+'_dir-pa_epi.nii.gz') )
    if not os.path.isfile(pa_file):
        pa_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-fmri_dir-pa_epi.nii.gz') )
        
        
    appa_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'fmap_appa.nii.gz'))
    topup_prefix = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_results_appa'))
    topup_field = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_field_appa'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'fmap_appa-dc.nii.gz'))
    
    
    # run topup
    if not os.path.isfile(topup_field + '.nii.gz'):
        
        utils.prettyOut("estimating field distortions for "+ task)
        
        # concatenate fmaps
        if not os.path.isfile(appa_file):
            utils.Cmd("fslmerge -t {} {} {}".format(appa_file, ap_file, pa_file)).run()

        # estimate distortion
        utils.Cmd("topup --verbose --imain=" + appa_file + " --datain=" + subj.study.paramsPath + "/epi_params.txt --config=" + subj.study.paramsPath+ "/b02b0.cnf --out=" + topup_prefix + " --fout=" + topup_field + " --iout=" + o_file).run()

        
        
    # clean up
    if os.path.isfile(o_file):
        utils.purge_paths(appa_file, o_file)





        
def apply_distortion_correction(subj, task, i_file, o_file, param_indx=1):
    """
    apply distortion correction to correct distortions for epi task    
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    i_file : TYPE
        DESCRIPTION.
    o_file : TYPE
        DESCRIPTION.
    param_indx : TYPE, optional
        DESCRIPTION. The default is 1.
        1,2 = AP; 3,4 = PA (in current params file)
        the index here is used to specify the acquisition encoding of the input and match it to the epi_params file.
        For our current setup, it will be '1' to specify that the entire thing was acquired A -> P
        
    Returns
    -------
    None.

    """
    utils.prettyOut("Applying distortion correction")
    if not os.path.isfile(o_file):
        
        topup_field = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_field_appa'))
        
        if not os.path.isfile(o_file):
            utils.Cmd("applytopup " + 
                      " --verbose " +
                      " --imain=" + i_file +
                      " --inindex=" + str(param_indx) + 
                      " --topup=" + topup_field +
                      " --datain=" + subj.study.paramsPath  + "/epi_params.txt"
                      " --method=jac " +
                      " --interp=spline" + 
                      " --datatype=short" + 
                      " --out=" + o_file).run()
    
            # remove negative values around the brain as a result of jac interpolation from dewarped output
            # make sure datatype INT
            utils.Cmd("fslmaths {0} -abs {0} -odt short".format(o_file)).run()        



    
    
def apply_distortion_correction_sbref(subj, task):
    
    sbref_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz'))
    apply_distortion_correction(subj, task, i_file=sbref_file, o_file=o_file)





def apply_distortion_correction_concat(subj, task):
    
    i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-uncorrected.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))
    apply_distortion_correction(subj, task, i_file=i_file, o_file=o_file)




    
def generate_task_mask(subj, task):
    
    """
    create mask on base(SBREF) image using dist corrected file

    """

    try:
        sbref_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz'))
        mask_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
        
        if not os.path.isfile(mask_file):
            utils.prettyOut(task + " : Make mask")
            utils.Cmd("3dAutomask -dilate 2 -erode 1 -prefix " + mask_file + " " + sbref_file).run()
        
    except:
        raise
        
    # apply this mask to sbref
    apply_task_mask(subj, task, sbref_file, os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz')))
    
    
    
    
    
def apply_task_mask(subj, task, i_file, o_file):
    """ 
    Apply task mask to i_file naming it o_file
    """
    if not os.path.isfile(i_file):
        raise ValueError("Input file not found")
        
    if not os.path.isfile(o_file):
        mask_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
        utils.Cmd( "fslmaths {} -mas {} {}".format(i_file, mask_file, o_file)).run()
            
            
            
            
    
def align_to_T1_BBR(subj, task, template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Boundary-Based Registration works better than ANTs to MNI template
    
    """
    
    anat_path = subj.derivatives_path(bidsType='anat')
    sbref_masked_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
    t1_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    t1_brain_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    t1_wmseg_path = os.path.join(anat_path, subj.bids_file_join('spm-wm_prob.nii.gz'))
    tx_for_pref = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-T1.aff'))
    tx_inv_pref = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'from-T1.aff'))
    
    


    try:
        
        if not os.path.isfile(tx_for_pref + '.mat'):   
            
            try:
                # do a standard flirt pre-alignment
                utils.Cmd("flirt -ref {} -in {} -dof 12 -omat {}".format(t1_brain_path, sbref_masked_file, os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')))).run()
        
                # now run the bbr
                utils.Cmd("flirt -ref {} -in {} -dof 12 -cost bbr -wmseg {} -init {} -omat {} -schedule {}/etc/flirtsch/bbr.sch".format(
                    t1_brain_path, sbref_masked_file, t1_wmseg_path,
                    os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')),
                    tx_for_pref + '.mat',
                    utils.Env.fsl_path
                    )).run()
            
            except:
                raise
                
            finally:
                # remove init mat
                utils.purge_paths(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')))
        
    
    
        # also save out inverse transform
        if not os.path.isfile(tx_inv_pref + '.mat'):
            utils.Cmd("convert_xfm -omat {} -inverse {}".format(tx_inv_pref + '.mat', tx_for_pref + '.mat')).run()
            
        
        # convert FSL transform to ITK format for ANTs
        if not os.path.isfile(tx_for_pref + '.mat.itk.txt'):
            utils.Cmd(utils.Env.c3d_path + "/c3d_affine_tool -ref {} -src {} {} -fsl2ras -oitk {}".format(t1_path, sbref_masked_file, tx_for_pref + '.mat', tx_for_pref + '.mat.itk.txt' )).run()
          
        
        # collapse the transformations to a final displacement field
        composite_warp_path = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-'+template_key+'_Composite.warp.nii.gz'))
        composite_in_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))
        if not os.path.isfile(composite_warp_path):
            utils.Cmd("antsApplyTransforms -d 3 -o '[{},1]' -t {} -t {} -r {}".format(composite_warp_path, composite_in_path, tx_for_pref + '.mat.itk.txt', utils.Env.template_lib[template_key]['brain'])).run()

    
    except:
        raise
        
    finally:
        utils.purge_paths(tx_for_pref + '_fast_*', tx_for_pref + '_init*',tx_for_pref + '.nii.gz')
    
    
    



def align_to_T1_AFNI(subj, task, cost_func="lpc", move='giant', cmass="-cmass cmass", opts="", purge=False):

    """
    Performs affine registration and outputs transform dwi2tlrc.aff12.1D

    Expects masked b0 image will already be available.

    :param subj: utils.Subj
    :param study: utils.Study
    :param cost_func: valid cost functions for 3dAllineate
    :param move:
    :param cmass:
    :param opts: Any other valid 3dAllineate flags to pass along
    :param purge: By default this function will not overwrite transforms as those may have been checked. Use this if you don't want to manually delete the previous transform.

    """
    anat_path = subj.derivatives_path(bidsType='anat')
    func_path = subj.derivatives_path(bidsType='func')
    tx_o_path = os.path.join(func_path, subj.bids_file_join('task-' + task.lower() + '-to-T1.aff12.1D'))
    sbref_masked_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
    
    
    # if alignment matrix exists and not told to blow it away
    if purge or not os.path.isfile(tx_o_path):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+task.lower()+" : Spatial Normalization" )

        # just so user understands
        if purge:
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " :Purging previous transform matrix")
            utils.purge_paths(tx_o_path)    


        #    create alignment workspace
        workspace_path = os.path.join(func_path, "regTemp")
        utils.make_path(workspace_path)

        # check for inputs
        if not os.path.isfile(sbref_masked_file):
            utils.mask_data(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz')),
                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz')),
                            sbref_masked_file)


        #    go in there cuz afni script is too dumb to use full paths correctly
        current_path = os.getcwd()
        os.chdir(workspace_path)

        #    do alignment, sadly afni script not made callable from python
        if (move) and (move in ['big', 'giant', 'ginormous']):
             move_str = "-"+move+"_move"
        else:
             move_str = ""

        utils.Cmd("align_epi_anat.py \
                -epi2anat \
                -anat " + os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz')) + " \
                -epi " + sbref_masked_file + " \
                -epi_base 0 \
                -epi_strip None \
                -volreg off \
                -deoblique on \
                -anat_has_skull no \
                -tshift off \
                -align_centers yes \
                "+ move_str +" \
                -cost "+ cost_func +" \
                -Allineate_opts '-source_automask+4'").run()

        
        import shutil        
        shutil.copy( subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked_al_mat.aff12.1D'), tx_o_path)

        # clean up
        os.chdir(current_path)
        shutil.rmtree(workspace_path, ignore_errors=True)


    return tx_o_path






def spatial_normalization(subj, task, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    
    Creates a final Nonlinear tramsform from functional to template.
    
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    try:
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+task.lower()+" : Spatial Normalization" )
        
        
        if space == 'TLRC':
            
            # make sure there is an affine transform
            # if target space is TLRC, use AFNI method for alignment to T1
            # if MNI, use BBR method            
            affine_tx = align_to_T1_AFNI(subj, task)
            
            anat_path = subj.derivatives_path(bidsType='anat')        
            o_path = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz'))
    
            if not os.path.isfile(o_path):
                utils.Cmd("3dNwarpCat -warp1 {} -warp2 {} -warp3 {} -prefix {} ".format(os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.WARP.nii.gz')),
                                                                                        os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.aff12.1D')),
                                                                                        affine_tx,
                                                                                        o_path)).run()
                
        else:
            affine_tx = align_to_T1_BBR(subj, task, template_key=template_key)
            

        
    except:
        raise








def apply_spatial_normalization(subj, task, i_file, o_path=None, o_file=None, master_path=None, interp='wsinc5', regrid=None, fwhm=None, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Spatially normalize i_file into o_path with option flags
    o_file will be i_file with appended flags options
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """


    
    
    
    try:

        # TODO: clean up this method, add other nice options like TLRC
        
        
        
        ####################################
        # MNI SPACE
        ####################################
        if space == 'MNI':  
            
            if not os.path.isfile(i_file):
                raise ValueError("Input file does not exist")
                
            # split input into path and file
            (i_path, i_file_name) = os.path.split(i_file)
    
            # reduce input file to root
            i_file_root = i_file_name.replace('.nii', '').replace('.gz','')        
            
            # use i_path as o_path if user did not specify
            if o_path == None:
                if o_file == None:
                    o_path = i_path
                    utils.make_path(o_path)
                else:
                    (o_path, o_file_name) = os.path.split(o_file)
                    utils.make_path(o_path)
                    
                    
            # build output filename
            if o_file == None:
                o_file_name = i_file_root
                if regrid != None:
                    o_file_name = o_file_name + '_{0}x{0}x{0}'.format(regrid)
        
                if fwhm != None:
                    o_file_name = o_file_name + '_GB{0}'.format(fwhm)
                    
                o_file_name = o_file_name + '_'+space+'.nii.gz'
            
            else:
                (o_path, o_file_name) = os.path.split(o_file)
                
            
            if master_path == None:
                master_path = utils.Env.template_lib[template_key]['brain']                 


    
            utils.Cmd("antsApplyTransforms -d 3 -o {} -t {} -r {} -i {} --output-data-type short --interpolation LanczosWindowedSinc --input-image-type 0".format(
                os.path.join(o_path, o_file_name),
                os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-'+template_key+'_Composite.warp.nii.gz')),
                master_path,
                os.path.join(i_path, i_file_name)
                )).run()
            
            # update NIFTI header
            from mayerlab.preprocessing.anat import afni
            afni.set_header_space(space='MNI', f_path=os.path.join(o_path, o_file_name))
    
            
    


        ####################################
        # TLRC SPACE
        ####################################
        if space == 'TLRC':
            
            if not os.path.isfile(i_file):
                raise ValueError("Input file does not exist")
                
            # split input into path and file
            (i_path, i_file_name) = os.path.split(i_file)
    
            # reduce input file to root
            i_file_root = i_file_name.replace('.nii', '').replace('.gz','')        
            
            # use i_path as o_path if user did not specify
            if o_path == None:
                if o_file == None:
                    o_path = i_path
                    utils.make_path(o_path)
                else:
                    (o_path, o_file_name) = os.path.split(o_file)
                    utils.make_path(o_path)
                
                
            # build output filename
            if o_file == None:
                o_file_name = i_file_root
                if regrid != None:
                    o_file_name = o_file_name + '_{0}x{0}x{0}'.format(regrid)
        
                if fwhm != None:
                    o_file_name = o_file_name + '_GB{0}'.format(fwhm)
                    
                o_file_name = o_file_name + '_'+space+'.nii.gz'
            
            else:
                (o_path, o_file_name) = os.path.split(o_file)
                
            
            
            
            anat_path = subj.derivatives_path(bidsType='anat')
            
            if master_path == None:
                if space == 'MNI':
                    master_path = os.path.join(anat_path, subj.bids_file_join('T1w_MNI.nii.gz'))
                else:
                    master_path = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.nii.gz'))
            
            
            
            if fwhm == None:
                if regrid == None:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            interp,
                                                                                                            os.path.join(o_path, o_file_name)
                                                                                                            )).run()
                    
                else:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -newgrid {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            str(regrid),
                                                                                                            interp,
                                                                                                            os.path.join(o_path, o_file_name)
                                                                                                            )).run()
                    
                
            else:
                
                if regrid == None:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            interp,
                                                                                                            os.path.join(o_path, 'tmp-'+o_file_name)
                                                                                                            )).run()                
                else:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -newgrid {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            str(regrid),
                                                                                                            interp,
                                                                                                            os.path.join(o_path, 'tmp-'+o_file_name)
                                                                                                            )).run()
    
                utils.Cmd("3dmerge -1blur_fwhm {} -doall -prefix {} {}".format(fwhm, os.path.join(o_path, o_file), os.path.join(o_path, 'tmp-'+o_file_name) )).run()
                utils.purge_paths(os.path.join(o_path, 'tmp-'+o_file_name))
                
            print("\nInput: " + os.path.join(i_path, i_file_name))    
            print("Output: " + os.path.join(o_path, o_file_name))
                
    except:
        raise








def apply_blur(subj, i_file, fwhm=None, o_path=None):
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

        
        if not os.path.isfile(i_file):
            raise ValueError("Input file does not exist")
            
        # split input into path and file
        (i_path, i_file) = os.path.split(i_file)

        # reduce input file to root
        i_file_root = i_file.replace('.nii', '').replace('.gz','')        
        
        # use i_path as o_path if user did not specify
        if o_path == None:
            o_path = i_path
            utils.make_path(o_path)
            
            
        # build output filename
        o_file = i_file_root
        if fwhm != None:
            o_file = o_file + '_GB{0}'.format(fwhm)
            
        o_file = o_file + '.nii.gz'
        
    
        utils.Cmd("3dmerge -1blur_fwhm {} -doall -prefix {} {}".format(fwhm, os.path.join(o_path, o_file), os.path.join(i_path, i_file) )).echo()
            
        print("\nOutput:" + os.path.join(o_path, o_file))
            
    except:
        raise



def clean_proc(subj, task):
    """
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_sbref-dc.nii.gz')),
    
    try:
        
        if os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
            
            utils.purge_paths( os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_dspk.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_tshft.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_2dreg*')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreg.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreglog.raw')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_log.mc')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_sbref.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_slice_timing.1D')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_preproc-uncorrected.nii.gz'))
                              )
                              
        
    except:
        raise







    

          
            
    
class qc:
    
    def run_mriqc(subj):
        """
        Run MRIQC for each task, each run as needed
        """
        
        proc_map = subj.get_subj_task_map()

        for task, run_number in proc_map.items():
            for run in range(1, run_number + 1):

                if not os.path.isfile(os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func', subj.bids_file_join('task-'+task.lower()+'_run-'+str(run)+'_bold.json'))):
        
                    # create workspace path outside of mriqc for correction group permissions
                    os.umask(0o007)
                    utils.make_path(os.path.join(subj.study.mriqcPath, 'workspace', 'func'))
                    
                    # must have this env variable set to have access to templates; does not work to include as --env 
                    os.environ['SINGULARITYENV_TEMPLATEFLOW_HOME'] = "/home/bidsapp/.cache/templateflow"
                    utils.Cmd("singularity run --cleanenv --containall --disable-cache --writable-tmpfs --bind /export:/export /export/research/analysis/human/amayer/shared/apps/containers/runtime/mayerlab_mriqc_0.16.1.sif "
                              + " --participant_label " + subj.ursi.short 
                              + " --session-id " + str(subj.visit)
                              + " --modalities bold"
                              + " --task-id " + task.lower()
                              + " --run-id " + str(run)
                              + " --work-dir " + subj.study.mriqcPath + "/workspace/func/" + subj.bids_file_join("task-"+task.lower()+"_run-"+str(run)+"_bold")
                              + " --no-sub "
                              + " --verbose "
                              + " --nprocs 6 "
                              + " --verbose-reports "
                              + " --float32 "
                              + " --ants-nthreads 6 "
                              + " --omp-nthreads 6 "
                              + subj.study.bidsPath + " "
                              + subj.study.mriqcPath + " participant "
                              ).run()
                    
                    # have to clean up after poldrack
                    utils.purge_paths(subj.study.mriqcPath + "/workspace/func/" + subj.bids_file_join("task-"+task.lower()+"_run-"+str(run)+"_bold"))

            try:
                # fix permissions on container task output for group access
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id) )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit) )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func') )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func', '*.json') )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_file_join('task-'+task.lower() + '*.html')) )

            except:
                pass




    def snapshot_mask(subj, task):
        
        
        """
        
        """
            
        if not os.path.isfile(os.path.join(subj.derivatives_path('func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask.axi.png'))):
            
            ulay_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, "preproc-dc.nii.gz'[0]'"))
            olay_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC Mask Snapshots" )
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+olay_file+" \
                    -prefix "+ os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask')) +" \
                    -set_xhairs OFF \
                    -montx 10 \
                    -monty 3 \
                    -label_string "+ subj.bids_file_join('task-' + task.lower(), 'mask') +" \
                    -pbar_posonly \
                    -cbar_ncolors 1").run()
                    
                    
            for view in ['axi', 'cor', 'sag']:
                utils.Cmd("convert " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask.' +view+'.png')) +" -gravity North -background YellowGreen -splice 0x18 -pointsize 18 -annotate +0+2 '"+task+" mask "+view+"'  " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask.' +view+'.png')) ).run()
                    
         
                
                
                
    def snapshot_spatial_normalization(subj, task, space='MNI'):
        
        
        """
        
        """
        if not os.path.isfile(os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space+'.jpg'))):

            src_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
            spa_norm_lay = os.path.join(subj.derivatives_path('func'), 'qc', 'sbref.masked.'+space+'.nii.gz')
            apply_spatial_normalization(subj, 
                                        task, 
                                        src_file, 
                                        o_file=spa_norm_lay,
                                        space=space)
    
            edge_first_image = os.path.join(subj.derivatives_path(bidsType='func'), 'qc', 'sbref.masked.'+space+'.edge.nii.gz')
            utils.Cmd("3dedge3 -overwrite  -input "+ spa_norm_lay +" -prefix " + edge_first_image).run()
            
            if space == 'MNI':
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/AFNI_MNI152_T1_2009c+tlrc.HEAD"
            else:
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/TT_N27+tlrc.HEAD"
                
                
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC '+space+' Snapshots" )
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+edge_first_image+" \
                    -prefix "+ os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space)) +" \
                    -label_string "+ subj.bids_file_join('task-' + task.lower(), space) + " \
                    -set_xhairs OFF \
                    -delta_slices 4 4 4 \
                    -montx 6 \
                    -monty 3 \
                    -opacity 9 \
                    -pbar_posonly \
                    -cbar_ncolors 6 \
                    -cbar_topval '' \
                    -cbar '1000=blue \
                            800=cyan \
                            600=rbgyr20_10 \
                            400=rbgyr20_08 \
                            200=rbgyr20_05 \
                            100=none \
                            0=none'").run()
    
    
            for view in ['axi', 'cor', 'sag']:
                utils.Cmd("convert " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) +" -gravity North -background YellowGreen -splice 0x18 -pointsize 18 -annotate +0+2 '"+task+" "+space+" "+view+"'  " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) ).run()
    
    
            #utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, spa_norm_lay, os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space)))).run()
            
            utils.purge_paths(spa_norm_lay,edge_first_image)

        
    def snapshot_spatial_normalization_old(subj, task, space='MNI'):
        
        
        """
        
        """
        if not os.path.isfile(os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space+'.jpg'))):

            src_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
            spa_norm_lay = os.path.join(subj.derivatives_path('func'), 'qc', 'sbref.masked.'+space+'.nii.gz')
            apply_spatial_normalization(subj, 
                                        task, 
                                        src_file, 
                                        o_file=spa_norm_lay,
                                        space=space)
    
            #edge_first_image = os.path.join(subj.derivatives_path(bidsType='func'), 'qc', 'sbref.masked.'+space+'.edge.nii.gz')
            #utils.Cmd("3dedge3 -overwrite  -input "+ spa_norm_lay +" -prefix " + edge_first_image).run()
            
            if space == 'MNI':
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/AFNI_MNI152_T1_2009c+tlrc.HEAD"
            else:
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/TT_N27+tlrc.HEAD"
                
                
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC '+space+' Snapshots" )
            
            #utils.Cmd("@chauffeur_afni \
                    #-ulay "+ulay_file+" \
                    #-olay "+edge_first_image+" \
                    #-prefix "+ os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space)) +" \
                    #-label_string "+ subj.bids_file_join('task-' + task.lower(), space) + " \
                    #-set_xhairs OFF \
                    #-delta_slices 4 4 4 \
                    #-montx 6 \
                    #-monty 3 \
                    #-opacity 9 \
                    #-pbar_posonly \
                    #-cbar_ncolors 6 \
                    #-cbar_topval '' \
                    #-cbar '1000=blue \
                            #800=cyan \
                            #600=rbgyr20_10 \
                            #400=rbgyr20_08 \
                            #200=rbgyr20_05 \
                            #100=none \
                            #0=none'").run()
    
    
            #for view in ['axi', 'cor', 'sag']:
                #utils.Cmd("convert " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) +" -gravity North -background YellowGreen -splice 0x18 -pointsize 18 -annotate +0+2 '"+task+" "+space+" "+view+"'  " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) ).run()
    
    
            utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, spa_norm_lay, os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space)))).run()
            
            utils.purge_paths(spa_norm_lay)
        
        
    def aggregate_qc_pdf(subj):
        
        
        """
        
        """    
        import glob
        
        
        func_path = subj.derivatives_path(bidsType='func')
        
        o_file = os.path.join(func_path, 'qc', subj.bids_file_join('func-QC.pdf'))
        if not os.path.isfile(o_file):
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC PDF" )
        
            os.chdir(os.path.join(func_path, 'qc'))
            
            files_for_pdf = ''
            jpg_list = glob.glob('*.jpg')
            if len(jpg_list):
                files_for_pdf = files_for_pdf + ' *.jpg'
            png_list = glob.glob('*.png')
            if len(png_list):
                files_for_pdf = files_for_pdf + ' *.png'
            if len(files_for_pdf):
                utils.Cmd("convert "+files_for_pdf+" -quality 100 " + o_file).run()
            
            
            
            
