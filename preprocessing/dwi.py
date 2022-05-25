"""@package dwi
Preprocessing module for diffusion imaging
Author: jling@mrn.org
DWI,DTI,DKI,etc.

Along with utils this script is mostly designed to work within a container.
"""

from mayerlab.preprocessing import utils
import os,errno




class pipelines:
    
    def __init__(self, pipeline_name):
        self.pipeline_name = pipeline_name
        

    def run(self, subj):
        """
        Intended for call as a container. If not in container call functions directly.
        """
        import inspect
        
        bFound_pipeline = False
        for f in inspect.getmembers(pipelines, predicate=inspect.isfunction):
            if f[0] == self.pipeline_name:
                bFound_pipeline = True
                f[1](self, subj)
        
        if not bFound_pipeline:
            raise ValueError('Unknown pipeline requested:' + self.pipeline_name)
            
            
            
            
    def standard(self, subj):
        """
        Defines standard preprocessing pipeline steps
    
        Parameters
        ----------
        subj
        
        Returns
        -------
        None.
    
        """
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+ self.pipeline_name +" Pipeline" + " : DWI Processsing")
        
        # avoid running on workstations that have one GPU as that is probably not appropriate
        # many workstations have a Quadro 600 that won't for this task
        if utils.Gpu().get_gpu_count() < 2:
            raise ValueError("This script requires a GPU server.")

        
        ###  concatenate runs
        concatenate_run_data(subj)

        ###  denoise data
        denoise.MPPCA(subj)
        
        ###  run EDDY for motion and susceptibility
        susceptibility_eddy_correction(subj, s2v=True, mbs=True, ol=False)

        ###  write out movement parameters as estimated with Eddy
        get_eddy_movement_data(subj)
        
        ###  generate qc metrics of EDDY output
        qc.run_eddy_quad(subj)
        
        ### 
        spatial_normalization(subj, space='MNI', template_key='mni_icbm152_nlin_asym_2009c')
        
        ### DTI
        model.generate_nonlinear_DT(subj, use_bvals='all')
        
        ### compile all QC images into PDF
        qc.aggregate_qc_pdf(subj, space='MNI')
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI Processsing Complete")
                
                
                
            
            
            
    def laptop(self, subj):
        """
        Defines lapto compatible preprocessing pipeline steps
    
        Parameters
        ----------
        subj
        
        Returns
        -------
        None.
    
        """

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+self.pipeline_name+" Pipeline" + " : DWI Processsing")
        
        
        # avoid running on workstations that have one GPU as that is probably not appropriate
        if utils.Gpu().get_gpu_count() < 2:
            raise ValueError("This script requires a GPU server.")
    
    
        # concatenate runs
        concatenate_run_data(subj)
    
        # run EDDY for motion and susceptibility
        susceptibility_eddy_correction(subj, s2v=True, mbs=True, ol=False)
    
        # write out movement parameters as estimated with Eddy
        get_eddy_movement_data(subj)
        
        # generate qc metrics
        qc.run_eddy_quad(subj)
        
        #
        spatial_normalization(subj, space='TLRC')

        # DTI
        model.generate_nonlinear_DT(subj, use_bvals='all')
        
        qc.snapshot_spatial_normalization(subj, space='TLRC')
        qc.aggregate_qc_pdf(subj, space='TLRC')

        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI Processsing Complete")


            







# def verify_phase_encoding(subj, fatal=False):

    # """

    # :param subj: utils.Subj

    # """

    # import json
    # import glob
    # import numpy as np
    
    
    # dwi_derv_path = subj.derivatives_path(bidsType='dwi')

    # #    create home for dwi output
    # try:
    #     os.makedirs(dwi_derv_path,0o770)
    # except OSError as exception:
    #     if exception.errno != errno.EEXIST:
    #         raise
    
    # dwi_raw_path = subj.bids_modality_path(bidsType='dwi')
    # if not os.path.isdir(dwi_raw_path):
    #     raise ValueError("Missing input:\n"+dwi_raw_path)
            
            
            
            
    # # concate dwi image files
    # o_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi_concat.nii.gz"))
    # if not os.path.isfile(o_file):
    #     utils.prettyOut("DWI : Concatenate Run Data")
    #     utils.Cmd( "fslmerge -t " +  o_file + ' ' + " ".join(sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.nii.gz'), key = lambda x: x.split('_')[-2])) ).run()
        
        
        
        
    # run_settings     = study.get_run_type_settings("DTI")

    # utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Phase Encoding Verification" )
    
    # for run_num in range(1,int(run_settings['run_count'])+1):
        
    #     print("DWI"+str(run_num))
        
    #     i_path = os.path.join(subj.path_append(study, proc_dir),"dwi"+str(run_num)+".json")
        
    #     if os.path.isfile(i_path):
    #         with open(i_path) as sc_fh:
    #             sc = json.load(sc_fh)
                
    #             if( (sc['ProtocolName'].split('_')[3] == 'ap') & (sc['PhaseEncodingDirection'] != 'j-') ):
    #                 print("PhaseEncodingDirection Error: A>P encoding not implemented.")
    #                 if fatal:
    #                     raise utils.CmdProcError("PhaseEncodingDirection Error")
    #             elif ( (sc['ProtocolName'].split('_')[3] == 'pa') & (sc['PhaseEncodingDirection'] != 'j') ):
    #                 print("PhaseEncodingDirection Error: P>A encoding not implemented.")
    #                 if fatal:
    #                     raise utils.CmdProcError("PhaseEncodingDirection Error")                    
    #             else:
    #                 print("PhaseEncodingDirection Correct")
    #     else:
    #         print("No BIDS sidecar file.")
    
    
    
    

def concatenate_run_data(subj):

    """
    Concatenates DWI data

    :param subj: utils.Subj

    """

    import glob
    import numpy as np
    
    
    dwi_derv_path = subj.derivatives_path(bidsType='dwi')

    #    create home for dwi output
    try:
        os.makedirs(dwi_derv_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    dwi_raw_path = subj.bids_modality_path(bidsType='dwi')
    if not os.path.isdir(dwi_raw_path):
        raise ValueError("Missing input:\n"+dwi_raw_path)
            
            
            
            
    # concate dwi image files
    o_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi-concat.nii.gz"))
    if not os.path.isfile(o_file):
        utils.prettyOut("DWI : Concatenate Run Data")
        utils.Cmd( "fslmerge -t " +  o_file + ' ' + " ".join(sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.nii.gz'), key = lambda x: x.split('_')[-2])) ).run()



    # concate bval files
    o_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi-concat.bvals"))
    if not os.path.isfile(o_file):
        utils.prettyOut("DWI : Concatenate Bval Data")        
        out_bvals_fh= open(o_file, "w")
        for i_bval_file in sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.bval'), key = lambda x: x.split('_')[-2]):
            with open(i_bval_file) as f:
                out_bvals_fh.write(f.readline().replace('\n', ' '))
                
        out_bvals_fh.write('\n')
        out_bvals_fh.close()

        
        
    # concate bvec files
    o_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi-concat.bvecs"))
    if not os.path.isfile(o_file):
        utils.prettyOut("DWI : Concatenate Bvec Data")        
        
        new = np.array([])
        for i_bvec_file in sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.bvec'), key = lambda x: x.split('_')[-2]):
            if new.shape[0] == 0:
                new = np.loadtxt(i_bvec_file)
            else:
                new = np.concatenate( (new, np.loadtxt(i_bvec_file) ), axis=1 )
                
        np.savetxt(o_file, new, delimiter=' ', fmt='%1.7f')

        
    # create rounded bval file
    round_bvals_file(subj)





def susceptibility_eddy_correction(subj, s2v=True, mbs=True, ol=False):

    """
    This corrects susceptibility-induced distortions with FSL EDDY
    using reversed phase-encoding polarities
    
    :param subj: utils.Subj
    :param study: utils.Study
    :param s2v: EDDY includes intra-volume (slice-to-vol) movement correction"
    :param mbs: EDDY includes susceptibility-by-movement correction
    :param ol: EDDY includes slice outlier replacement
    """
    
    dwi_path = subj.derivatives_path(bidsType='dwi')
    edc_o_path = os.path.join(dwi_path, subj.bids_file_join("dwi-edc.nii.gz"))
    if not os.path.isfile(edc_o_path):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Eddy current and motion correction" )
        
        # verify params input first
        file_path_verify = [os.path.join(subj.study.paramsPath, 'dwi_params.txt'), os.path.join(subj.study.paramsPath, 'dwi_index.txt')]
        for fp in file_path_verify:
            if not os.path.isfile(fp):
                raise utils.CmdProcError("Missing Param input: " + fp)


        # get b0 volumes for topup
        b0_filepath = os.path.join(dwi_path, subj.bids_file_join("dwi-topup_b0.nii.gz"))
        fmap_to_topup_b0(subj, b0_filepath)          
        
        # run topup
        run_topup(subj, b0_filepath, dwi_path)
        
        # create dwi mask and toput mean image
        create_dwi_mask(subj)
        # qc snapshots of topup
        qc.snapshot_b0(subj, b0_filepath)        
        # run EDDY
        run_eddy_cuda(subj, s2v=s2v, mbs=mbs, ol=ol)







def generate_fsl_eddy_slspec(i_path_json, o_path):
    """
    Writes spec file using dcm2niix json file.
    Currently needed by eddy quad for s2v qc
    
    i_path_json: path to the JSON file that will be used as source
    o_path: full path and filename of output
    """
    import json
    import numpy as np
    
    # read json file
    with open(i_path_json) as f:    
        data = json.load(f)
        # get slice timing value
        slice_timing = np.array(data['SliceTiming'], dtype=float)
        # get mb factor
        mb = np.array(data['MultibandAccelerationFactor'], dtype=int)
        
    
    # sort to get index order
    sindx = np.argsort(slice_timing, axis=0, kind='mergesort')
    
    # generate slspec matrix
    slspec = np.reshape(sindx, (mb, int(len(sindx)/mb)), order='F')

    # save to o_file
    np.savetxt(o_path, slspec.T, fmt='%.d')
    
    



def round_bvals_file(subj):

    """
    Rounds b-values so they can be used with packages that expect that.
    Function uses subject and study input to find current bvals file.
    From that it creates a file in the same place and name except ".rounded"
    is added to end of filename.


    :param subj: utils.Subj
    :param study: utils.Study
    :returns path to rounded bvals file
    
    See also: DIPY function

    """
    
    dwi_path = subj.derivatives_path(bidsType='dwi')
    in_bvals_path = os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals"))
    out_bvals_path = os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals-rounded"))
    
    if not os.path.isfile(out_bvals_path):
        if not os.path.isfile(in_bvals_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvals_path)

        out_bvals_fh= open(out_bvals_path, "w")
        with open(in_bvals_path) as f:
            for line in f:
                for bval in line.strip().split(' '):
                    bval_i = int(round(float(bval+".0")/100)*100)
                    out_bvals_fh.write(str(bval_i)+" ")

        out_bvals_fh.close()

    return out_bvals_path





def get_bvecs_by_bvals(subj, study, bvalue):

    """
    Returns space-delimited bvec rows given a b-value

    :param subj: utils.Subj
    :param study: utils.Study
    :param bvalue: restrict to a b-value

    """

    # we need rounded bvals file. Just call this function. It will check if it needs to create the file
    in_bvals_path = round_bvals_file(subj,study)

    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)


    in_bvecs_path = subj.path_append(study, 'dwi', subj.id_fn(study,"dwi.concat.bvecs.adj"))
    if not os.path.isfile(in_bvecs_path):
        in_bvecs_path = subj.path_append(study, 'dwi', subj.id_fn(study,"dwi.concat.bvecs"))
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)

    #cmd = utils.Cmd("cat " + in_bvals_path + " | tr ' ' '\n' | paste -d' ' - "+in_bvecs_path+" | grep '^"+bvalue+"' | cut -d' ' -f2-4")
    import subprocess

    cmd_txt = "cat " + in_bvals_path + " | tr ' ' '\\n' | paste -d' ' - " + in_bvecs_path + " | grep '^"+str(bvalue)+"' | cut -d' ' -f2-4"
    process = subprocess.Popen(cmd_txt,stdout=subprocess.PIPE, shell=True)
    cmd_result = process.stdout.read()

    return cmd_result





def run_eddy_cuda(subj, s2v=True, mbs=True, ol=False):
    """
    Runs FSL EDDY through Singularity, Docker, or locally if inside a container
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    study : TYPE
        DESCRIPTION.
    :param s2v: EDDY includes intra-volume (slice-to-vol) movement correction"
    :param mbs: EDDY includes susceptibility-by-movement correction
    :param ol: EDDY includes slice outlier replacement
    
    Returns
    -------
    None.

    """

    import shutil
    from math import ceil
    
    dwi_path = subj.derivatives_path(bidsType='dwi')
    edc_o_path = os.path.join(dwi_path, subj.bids_file_join("dwi-edc.eddy_command_txt"))
    if not os.path.isfile(edc_o_path):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Eddy ...Smoke, if ya got 'em" )

        try:
            
            imain=os.path.join(dwi_path, subj.bids_file_join("dwi-concat.nii.gz"))
            mask=os.path.join(dwi_path, subj.bids_file_join("dwi-edc_mask.nii.gz"))
            acqp=os.path.join(subj.study.paramsPath, 'dwi_params.txt')
            index=os.path.join(subj.study.paramsPath, 'dwi_index.txt')
            bvecs=os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvecs"))
            bvals=os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals-rounded"))
            topup=os.path.join(dwi_path, "topup_results")
            out=os.path.join(dwi_path, subj.bids_file_join("dwi-edc"))
            
            # get dwi json sidecar for slicetiming
            json=get_json_file_list(subj, indx=0)
            
            #  mporder: for a multi-band acquisition N would be the number of slices divided by the multi-band factor. 
            # If one wants to do slice-to-vol motion correction --mporder should be set to an integer value greater than 0 and less than the number of excitations in a volume. 
            # Only when --mporder > 0 will any of the parameters prefixed by --s2v_ be considered. The larger the value of --mporder, the more degrees of freedom for modelling movement. 
            # If --mporder is set to N-1, where N is the number of excitations in a volume, the location of each slice/MB-group is individually estimated. We don't recommend going that 
            # high and in our tests we have used values of N/4 -- N/2.
            # Jesper Andersson <jesper.andersson@NDCN.OX.AC.UK>
            dwi_param_keys = ["MultibandAccelerationFactor","SliceTiming"]
            dwi_param_dict = get_json_param_value_dict(json, dwi_param_keys)
            if len(dwi_param_dict.keys()) != len(dwi_param_keys):
                raise ValueError("get_json_param_value_dict() json key not found")
            
            mb = dwi_param_dict["MultibandAccelerationFactor"]
            slice_num = len(dwi_param_dict["SliceTiming"])
            mporder = str(ceil(slice_num/mb/4))

            # NVS 310 doesn't reveal status so we cannot choose a gpu
            use_gpu = utils.Gpu().get_available(55.0)
            if use_gpu != None:
                os.environ['SINGULARITYENV_CUDA_VISIBLE_DEVICES'] = use_gpu
    
    
            # if running inside of a container we will call eddy_cuda directly as it should be in the container
            # if running outside of a container, call the appropriate external eddy_cuda container
            if utils.Env.running_in_container():
                eddy_cmd_str = "eddy_cuda "
            else:
                if shutil.which("/usr/local/bin/singularity") == None:
                    eddy_cmd_str = "docker run -it --rm --runtime=nvidia -e NVIDIA_VISIBLE_DEVICES="+use_gpu+" --user `id -u`:`id -g` --mount type=bind,source=/export,target=/export mayerlab/eddy_cuda:0.0.6 "
                else:
                    eddy_cmd_str = "singularity run --nv --cleanenv --containall --disable-cache --bind /export:/export /export/research/analysis/human/amayer/shared/apps/containers/runtime/mayerlab_eddy_cuda_0.0.6.sif "
                
            
            # append basic run inputs
            # imain, mask, acqp, index, bvecs, bvals, topup, out, json)
            eddy_cmd_str = eddy_cmd_str + "--imain={} --mask={} --acqp={} --index={} --bvecs={} --bvals={} --topup={} --out={} --json={} --flm=quadratic --slm=none --fwhm=10,0,0,0,0 --niter=5 --interp=spline --resamp=jac --nvoxhp=1000 --ff=10 --cnr_maps --data_is_shelled --verbose ".format(
                imain, mask, acqp, index, bvecs, bvals, topup, out, json)
            
            
            
            #  hidden command sometimes useful
            #  --dont_mask_output=True
            
            # if requested, add other features to eddy
            if s2v:
                print("    +++ EDDY includes intra-volume (slice-to-vol) movement correction")
                eddy_cmd_str = eddy_cmd_str + "--mporder="+mporder+" --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear "
            if mbs:
                print("    +++ EDDY includes susceptibility-by-movement correction")
                eddy_cmd_str = eddy_cmd_str + "--estimate_move_by_susceptibility --mbs_niter=10 --mbs_lambda=10 --mbs_ksp=10 "
            if ol:
                print("    +++ EDDY includes slice outlier replacement")
                eddy_cmd_str = eddy_cmd_str + "--repol --ol_nstd=4 --ol_nvox=250 --ol_type=sw "
                
            
            # run eddy
            utils.Cmd(eddy_cmd_str).run()

            # remove negative values around the brain as a result of jac interpolation
            # convert float datatype that comes out of EDDY
            if os.path.isfile( os.path.join(dwi_path, subj.bids_file_join("dwi-edc.nii.gz")) ):
                utils.Cmd("fslmaths " + os.path.join(dwi_path, subj.bids_file_join("dwi-edc.nii.gz")) + " -abs " + os.path.join(dwi_path, subj.bids_file_join("dwi-edc")) + " -odt short").run()
            
            # fix group permissions on container output
            utils.set_lab_file_permissions(os.path.join(dwi_path, "*-edc.*"))

        except:
            # clean eddy remants on failure
            utils.purge_paths(os.path.join(dwi_path, subj.bids_file_join("dwi-edc.eddy_*")))
            raise






def get_eddy_movement_data(subj, out_units='radians'):

    """
    extracts and cleans dwi motion parameters output from FSL EDDY
    
    
    
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#HelpOnRotBvec
    Understanding eddy output:
        eddy_output.eddy_parameters : This is a text file with one row for each volume and one column for each parameter. 
                                      The first six columns correspond to subject movement starting with three translations followed by three rotations. 
                                      The remaining columns pertain to the EC-induced fields and the number and interpretation of them will depend of which EC model was specified.


        eddy_output.eddy_movement_rms : A summary of the "total movement" in each volume is created by calculating the displacement of each voxel and then averaging the squares of those displacements across all intracerebral voxels (as determined by --mask and finally taking the square root of that. The file has two columns where the first contains the RMS movement relative the first volume and the second column the RMS relative the previous volume.
    
    
    
    Put output in a format similar to 3dvolreg to avoid confusion.
    i.e.
    roll pitch yaw dS dL dP

                    roll  = rotation about the I-S axis }
                    pitch = rotation about the R-L axis } degrees CCW
                    yaw   = rotation about the A-P axis }
                      dS  = displacement in the Superior direction  }
                      dL  = displacement in the Left direction      } mm
                      dP  = displacement in the Posterior direction }

    ASSUMING RPI
        rot_z rot_x rot_y trans_z trans_x trans_y
        5     3     4     2       0       1


    R = Y X Z where Z is rotation about the I-S axis, X is rotation about the R-L axis, and Y is rotation about the A-P axis

    this does not give the same output layout as 3dvolreg but puts "edc" in file name
    Rotations in this file are given in radians unlike 3dvolreg which is in degrees

    :param subj: utils.Subj

    """

    dwi_path = subj.derivatives_path(bidsType='dwi')
    
    out_file_path = os.path.join(dwi_path, subj.bids_file_join("dwi-edc.Movement.Regressor.1D"))
    if not os.path.isfile(out_file_path):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Extract motion parameters" )
        
        #    if not file, complain
        in_file_path = os.path.join(dwi_path, subj.bids_file_join("dwi-edc.eddy_parameters"))
        if not os.path.isfile(in_file_path):
            raise utils.CmdProcError("Missing " + in_file_path)
    
        else:
            import numpy as np
    
            #    read rigid body alignment columns from eddy data file
            dwi_mov_data = np.loadtxt(in_file_path)[:,0:6]
            #    i don't like how reference image line is sometimes really small numbers so i force them to zero
            dwi_mov_data[0,:] = [0,0,0,0,0,0]
            
            if out_units == 'degrees':
                # convert radians to degrees
                dwi_mov_data[:,3:] = utils.Convert.Motion.radiansToDegrees(dwi_mov_data[:,3:])
                
            #    reorder; see notes above
            #    write out; place "edc" in file name to denote that this is not 3dvolreg output. e.g. rotations are in radians
            
            np.savetxt(out_file_path, dwi_mov_data[:,(5,3,4,2,0,1)], delimiter=' ', fmt='%1.6f')
    





        
        
def align_to_T1_MNI(subj, template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Boundary-Based Registration works better than ANTs to MNI template
    
    """
    
    anat_path = subj.derivatives_path(bidsType='anat')
    t1_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    t1_brain_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    t1_wmseg_path = os.path.join(anat_path, subj.bids_file_join('spm-wm_prob.nii.gz'))
    
    dwi_path = subj.derivatives_path(bidsType='dwi')
    b0_input_file = os.path.join(dwi_path, 'topup.dewarped_b0.Tmean.brain.nii.gz')
                                        
    tx_for_pref = os.path.join(dwi_path, subj.bids_file_join('dwi-to-T1.aff'))
    tx_inv_pref = os.path.join(dwi_path, subj.bids_file_join('dwi-from-T1.aff'))
    

    try:
        
        if not os.path.isfile(tx_for_pref + '.mat'):   
            
            try:
                
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Spatial Normalization")
                
                # make sure we haved masked input
                utils.mask_data(os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz"),
                                os.path.join(dwi_path, subj.bids_file_join("dwi-edc_mask.nii.gz")),
                                b0_input_file)
                
                # do a standard flirt pre-alignment
                print("   +++ initial transform")
                utils.Cmd("flirt -ref {} -in {} -dof 12 -omat {}".format(t1_brain_path, b0_input_file, os.path.join(dwi_path, subj.bids_file_join('dwi', 'init.mat')))).run()
        
                # now run the bbr
                print("   +++ BBR transform")
                utils.Cmd("flirt -ref {} -in {} -dof 12 -cost bbr -wmseg {} -init {} -omat {} -schedule {}/etc/flirtsch/bbr.sch".format(
                    t1_brain_path, b0_input_file, t1_wmseg_path,
                    os.path.join(dwi_path, subj.bids_file_join('dwi', 'init.mat')),
                    tx_for_pref + '.mat',
                    os.environ['FSLDIR']
                    )).run()
            
            except:
                raise
                
            finally:
                # remove init mat
                utils.purge_paths(os.path.join(dwi_path, subj.bids_file_join('dwi', 'init.mat')))
        
    
    
        # also save out inverse transform
        if not os.path.isfile(tx_inv_pref + '.mat'):
            print("   +++ inverse transform")
            utils.Cmd("convert_xfm -omat {} -inverse {}".format(tx_inv_pref + '.mat', tx_for_pref + '.mat')).run()
            
        
        # convert FSL transform to ITK format for ANTs
        if not os.path.isfile(tx_for_pref + '.mat.itk.txt'):
            print("   +++ convert transform for ANTs")
            utils.Cmd("c3d_affine_tool -ref {} -src {} {} -fsl2ras -oitk {}".format(t1_path, b0_input_file, tx_for_pref + '.mat', tx_for_pref + '.mat.itk.txt' )).run()
          
        
        # collapse the transformations to a final displacement field
        composite_warp_path = os.path.join(dwi_path, subj.bids_file_join('dwi-to-'+template_key+'_Composite.warp.nii.gz'))
        composite_in_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))
        if not os.path.isfile(composite_warp_path):
            print("   +++ composite transform")
            os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "22"
            utils.Cmd("antsApplyTransforms -d 3 -o '[{},1]' -t {} -t {} -r {}".format(composite_warp_path, composite_in_path, tx_for_pref + '.mat.itk.txt', utils.Env.template_lib[template_key]['brain'])).run()

    
    except:
        raise
        
    finally:
        utils.purge_paths(tx_for_pref + '_fast_*', tx_for_pref + '_init*',tx_for_pref + '.nii.gz')
    
    
    
    
        
        
        
def align_to_T1_TLRC(subj, cost_func="lpc", move='giant', cmass="-cmass cmass", opts="", purge=False):

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
    dwi_path = subj.derivatives_path(bidsType='dwi')
    tx_o_path = os.path.join(dwi_path, subj.bids_file_join("dwi-to-T1.aff12.1D"))
    
    
    # if alignment matrix exists and not told to blow it away
    if purge or not os.path.isfile(tx_o_path):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Spatial Normalization")
        
        import shutil

        # just so user understands
        if purge:
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " :Purging previous transform matrix")
            utils.purge_paths(os.path.join(dwi_path, subj.bids_file_join("dwi-to-T1.aff12.1D")))    


        #    create alignment workspace
        workspace_path = os.path.join(dwi_path, "regTemp")
        utils.make_path(workspace_path)

        # check for inputs
        b0_input_file = os.path.join(workspace_path, "topup.dewarped_b0.Tmean.sksp.nii.gz")
        utils.mask_data(os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz"),
                        os.path.join(dwi_path, subj.bids_file_join("dwi-edc_mask.nii.gz")),
                        b0_input_file)


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
                -epi " + b0_input_file + " \
                -epi_base 0 \
                -epi_strip None \
                -volreg off \
                -deoblique on \
                -anat_has_skull no \
                -tshift off \
                -cmass cmass "+ move_str +" \
                -cost "+ cost_func +" \
                -Allineate_opts '-source_automask+4'").run()


        shutil.copy("topup.dewarped_b0.Tmean.sksp_al_mat.aff12.1D", os.path.join(dwi_path, subj.bids_file_join("dwi-to-T1.aff12.1D")))

        # clean up
        os.chdir(current_path)
        shutil.rmtree(workspace_path, ignore_errors=True)






def spatial_normalization(subj, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    
    Parameters
    ----------
    subj : 
    space : Target space to use.
    template_key : reference name for template in that space to use as target
                   see utils.Env.template_lib for available templates
    """

    try:
        
        if space == 'MNI':
            
            # make sure there is an affine transform
            align_to_T1_MNI(subj, template_key=template_key)


        
        elif space == 'TLRC':
            # make sure there is an affine transform
            align_to_T1_TLRC(subj)
            
            anat_path = subj.derivatives_path(bidsType='anat')
            dwi_path = subj.derivatives_path(bidsType='dwi')
            
            o_path = os.path.join(dwi_path, subj.bids_file_join('dwi-to-'+space+'.NWARP.nii.gz'))
            if not os.path.isfile(o_path):
                utils.Cmd("3dNwarpCat -warp1 {} -warp2 {} -warp3 {} -prefix {} ".format(os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.WARP.nii.gz')),
                                                                                        os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.aff12.1D')),
                                                                                        os.path.join(dwi_path, subj.bids_file_join('dwi-to-T1.aff12.1D')),
                                                                                        o_path)).run()
        
        else:
            raise ValueError('spatial_normalization requires a valid space')
            
        
            
        # do qc snapshots
        qc.snapshot_spatial_normalization(subj, space=space, template_key=template_key)
        
        
        
    except:
        raise








def apply_spatial_normalization(subj, i_file, o_path=None, master_path=None, interp='wsinc5', regrid=None, fwhm=None, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    """

    try:

        if not os.path.isfile(o_path):
            
            dwi_path = subj.derivatives_path(bidsType='dwi')
            
            if space == 'MNI':

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
    
                os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "22"
                
                utils.Cmd("antsApplyTransforms -d 3 -o {} -t {} -r {} -i {} --output-data-type short --interpolation LanczosWindowedSinc --input-image-type 0".format(
                    o_path,
                    os.path.join(dwi_path, subj.bids_file_join('dwi-to-'+template_key+'_Composite.warp.nii.gz')),
                    master_path,
                    os.path.join(i_path, i_file_name)
                    )).run()
                
                
                # update NIFTI header
                from mayerlab.preprocessing.anat import afni
                afni.set_header_space(space='MNI', f_path=os.path.join(o_path))
                
                
                
                
            elif space == 'TLRC':
                
                if master_path == None:
                    #master_path = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.'+space+'.nii.gz'))
                    master_path = utils.Env.template_lib['TT_N27']['head']
                        
                utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(i_file,
                                                                                                       os.path.join(dwi_path, subj.bids_file_join('dwi-to-'+space+'.NWARP.nii.gz')),
                                                                                                       master_path,
                                                                                                       interp,
                                                                                                       o_path
                                                                                                       )).run()
            else:
                raise ValueError('spatial_normalization requires a valid space')
            
    except:
        raise





def index_by_pe(subj, bvalue=0, pe='all'):
    """
    Returns slicing index appropriate for extracting volumes by bvalue and or phase encoding
    """
    import numpy as np
    
    
    # read bvals and transpose
    dwi_path = subj.derivatives_path(bidsType='dwi')        
    bvals_path = os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals"))
    if not os.path.isfile(bvals_path):
        raise utils.CmdProcError("Missing bvec file for determining b0 volumes")
            
    # round bvalues because they will all be off a little 
    bvalues = np.loadtxt(bvals_path, dtype=int)
    bvalues = (np.round(bvalues/100.0)*100).astype(int)

    
    if pe == 'all':
        return list(np.where( bvalues == bvalue ).astype(str))
    else:
        # read dwi_index.txt
        dwi_index_path = os.path.join(subj.study.paramsPath, 'dwi_index.txt')
        phase_encoding = np.loadtxt(dwi_index_path, dtype=int)
        
        if pe.upper() == 'AP':
            return list(np.intersect1d(np.where(bvalues == bvalue), np.where(phase_encoding == 1), return_indices=False).astype(str))
        else:
            return list(np.intersect1d(np.where(bvalues == bvalue), np.where(phase_encoding == 2), return_indices=False).astype(str))
        

    


def fmap_to_topup_b0(subj, b0_filepath):
    """
    Looks for ap-pa dwi fmaps. 
    """
    
    dwi_raw_path = subj.bids_modality_path(bidsType='fmap')
    ap_file = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-ap_epi.nii.gz"))
    pa_file = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-pa_epi.nii.gz"))


    if not os.path.isfile(b0_filepath):    
        if os.path.isfile(ap_file) and os.path.isfile(pa_file):
            utils.Cmd("fslmerge -t {} {} {} ".format(b0_filepath, ap_file, pa_file) ).run()
        else:
            extract_dmap(subj)
            utils.Cmd("fslmerge -t {} {} {}".format(b0_filepath, ap_file, pa_file) ).run()
    
    
    
    
    
def extract_volume_by_index(in_file_path, out_file_path, b0_indices):

    """
    collects single AP PA pair for topup    
    """
    # check if already exists
    if not os.path.isfile( out_file_path ):
        utils.Cmd( "3dcalc -a " + in_file_path + "'[" + ",".join(b0_indices) + "]' -expr 'a' -prefix " + out_file_path ).run()
        





def extract_all_b0_to_volume(subj, out_file):

    """
    collects all the unweighted images from a raw dwi file and creates a new volume with these only

    :param subj: utils.Subj
    :param out_file Full file path for b0 output volume

    """

    dwi_path = subj.derivatives_path(bidsType='dwi')        
    if not os.path.isfile( out_file ):
        bvals_path = os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals"))
        if not os.path.isfile(bvals_path):
            raise utils.CmdProcError("Missing bvec file for determining b0 volumes")

        with open(bvals_path) as f:
            bvals_list = f.readline().strip().split()

        extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), 
                                out_file, 
                                [ str(i) for i in range(len(bvals_list)) if bvals_list[i] == '0' ])
        



def get_best_b0_index(subj, pe=None):
        """
        Use 3dTqual to select B0 volume
        """
        
        import subprocess
        
        if pe == None:
            raise utils.CmdProcError("Specify PE direction")
            
        # get all pe
        pe_indx = index_by_pe(subj, bvalue=0, pe=pe)
        if not len(pe_indx):
            raise utils.CmdProcError("Could not extract B=0")
        
        try:
            
            # work path
            dwi_path = subj.derivatives_path(bidsType='dwi')  
            # get all B0 for PE
            pe_b0_file = os.path.join(dwi_path, 'tmp_dwi-b0_'+pe+'.nii.gz')
            extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), pe_b0_file, pe_indx)
            
            # find b0 image which deviates least from the median
            # TODO:add cost function to make sure AP and PA are closest in space in addition to without motion?
            cmd_result = subprocess.run(["3dTqual", "-automask",  pe_b0_file], capture_output=True, text=True)
            qa_list =cmd_result.stdout.split()
            best_indx = qa_list.index(min(qa_list))
            
            return pe_indx[best_indx]    
        
        except:
            raise
            
        finally:
            utils.purge_paths( os.path.join( dwi_path, 'tmp_dwi-b0*') )
            
            
            
            
            
def extract_dmap(subj):

    """
    get best b=0 to volumes from dwi input into a dmap file
    This would be appropriate when fmaps are not available
    """
    
    
    dwi_raw_path = subj.bids_modality_path(bidsType='fmap')
    # expected dmap files
    ap_fmap_file = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-ap_epi.nii.gz"))
    pa_fmap_file = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-pa_epi.nii.gz"))
    
    dwi_path = subj.derivatives_path(bidsType='dwi')  
    
    
    # build separate AP and PA dmaps
    if not os.path.isfile(ap_fmap_file):
        print("\n   +++ extract AP fmap from concatenated volume B0s")
        try:
            utils.make_path(dwi_raw_path)
            # use best B0 as fmap
            extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), ap_fmap_file, [get_best_b0_index(subj, pe='AP')])
            # extracted b0 to fmap is going into rawdata which needs to comply with BIDS spec. Make that happy by creating JSON file. Could also just put the fmaps elsewhere but this is more consistent
            json_in_path = get_json_file_list(subj, substr='dir-AP')[0]
            json_out_path = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-ap_epi.json"))
            clone_json_file_params(json_in_path, json_out_path, params=['PhaseEncodingDirection','TotalReadoutTime'])

        except:
            raise
        
    if not os.path.isfile(pa_fmap_file): 
        print("\n   +++ extract PA fmap from concatenated volume B0s")
        try:
            utils.make_path(dwi_raw_path)
            # use best B0 as fmap
            extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), pa_fmap_file, [get_best_b0_index(subj, pe='PA')])
            # extracted b0 to fmap is going into rawdata which needs to comply with BIDS spec. Make that happy by creating JSON file. Could also just put the fmaps elsewhere but this is more consistent
            json_in_path = get_json_file_list(subj, substr='dir-PA')[0]
            json_out_path = os.path.join(dwi_raw_path, subj.bids_file_join("acq-dwi_dir-pa_epi.json"))
            clone_json_file_params(json_in_path, json_out_path, params=['PhaseEncodingDirection','TotalReadoutTime'])    
            
        except:
            raise
            
    # # write out to file which b0 indices we used. 
    # with open(os.path.join( dwi_path, subj.bids_file_join('dwi-topup_b0_list.txt') ), 'w') as f:
    #     f.write("AP volume index:"+ap_best_indx+"\n")
    #     f.write("PA volume index:"+pa_best_indx+"\n")

        






def extract_topup_b0_volume(subj, out_file):

    """
    get best b=0 to volumes from dwi input
    This would be appropriate when fmaps are not available
    """

    if not os.path.isfile(out_file):
            
        import subprocess
        
        # build separate AP and PA ouputs
        # get best AP and PA
        ap_indx = index_by_pe(subj, bvalue=0, pe='AP')
        if not len(ap_indx):
            raise utils.CmdProcError("Could not extract AP B=0")
        
        
        dwi_path = subj.derivatives_path(bidsType='dwi')  
           
        ap_file = os.path.join(dwi_path, 'tmp_dwi-b0_AP.nii.gz')
        extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), ap_file, ap_indx)
        ap_reg_file = os.path.join(dwi_path, 'tmp_dwi-b0_AP-reg.nii.gz')
        cmd_result = subprocess.run(["3dTqual", "-automask",  ap_reg_file], capture_output=True, text=True)
        qa_list =cmd_result.stdout.split()
        best_indx = qa_list.index(min(qa_list))
        ap_best_indx = ap_indx[best_indx]
            
        pa_indx = index_by_pe(subj, bvalue=0, pe='PA')
        if not len(ap_indx):
            raise utils.CmdProcError("Could not extract PA B=0")
            
        pa_file = os.path.join(dwi_path, 'tmp_dwi-b0_PA.nii.gz')
        extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')), pa_file, pa_indx)
        pa_reg_file = os.path.join(dwi_path, 'tmp_dwi-b0_PA-reg.nii.gz')
        cmd_result = subprocess.run(["3dTqual", "-automask",  pa_reg_file], capture_output=True, text=True)
        qa_list =cmd_result.stdout.split()
        best_indx = qa_list.index(min(qa_list))
        pa_best_indx = pa_indx[best_indx]
    
    
        utils.purge_paths( os.path.join( dwi_path, 'tmp_dwi-b0*') )
        
        with open(os.path.join( dwi_path, subj.bids_file_join('dwi-topup_b0_list.txt') ), 'w') as f:
            f.write("AP volume index:"+ap_best_indx+"\n")
            f.write("PA volume index:"+pa_best_indx+"\n")
            
        
        #   extract to topup_b0
        extract_volume_by_index(os.path.join(dwi_path, subj.bids_file_join('dwi-concat.nii.gz')),
                                out_file,
                                [ ap_best_indx, pa_best_indx ])
        






def run_topup(subj, b0_filepath, out_file_path):

    """ 
    run topup for susceptibility estimate
    
    In general I don’t recommend using more than a single pair of b0 images to run topup. Using more doesn’t seem to improve the estimation, 
    and only increases the execution time. One important reason that we still recommend acquiring more than one b0 or each PE-direction is 
    to increase ones changes to obtain at least one b0 volume that is not affected by “within volume subject movement”. The recommendation 
    is to manually inspect the b0 images, and then select the “best” for extraction for topup. Most of the time there will not be any problem 
    with subject movement in either of them, and then we recommend using the first. Jesper Andersson <jesper.andersson@NDCN.OX.AC.UK>
    
    """
    
    #    WARNING: Eddy cannot handle periods in the results file names
    if not os.path.isfile(os.path.join(out_file_path, "topup_results_fieldcoef.nii.gz") ):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Running Topup")

        # determine highest image multiple to speed up topup
        cnf_choice_str = get_highest_image_multiple(b0_filepath)
        # estimate distortion
        utils.Cmd("topup --verbose "
                + " --imain=" + b0_filepath
                + " --datain=" + os.path.join(subj.study.paramsPath, 'dwi_params.txt')
                + " --config=b02b0_"+cnf_choice_str+".cnf"
                + " --out=" + os.path.join(out_file_path, "topup_results")
                + " --iout=" + os.path.join(out_file_path, "topup.dewarped_b0")
                + " --logout=" + os.path.join(out_file_path, "topup_log") ).run()

        
        # remove negative values around the brain as a result of jac interpolation from dewarped output
        # also set datatype back to INT
        utils.Cmd("fslmaths {0} -abs {0} -odt short".format(os.path.join(out_file_path, "topup.dewarped_b0"))).run()





def create_dwi_mask(subj):

    """ 
    create dwi mask used for many pre and post commands
    TODO: fill in suscep holes??? check against edc
    
    """
    
    dwi_path = subj.derivatives_path(bidsType='dwi')  
    mask_file = os.path.join(dwi_path, subj.bids_file_join('dwi-edc_mask.nii.gz'))
        
        
    if not os.path.isfile(mask_file):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Create Mask")

        #    average susceptibility distortion corrected b0 volumes
        #    also set datatype back to INT
        if not os.path.isfile(os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz" )):
            utils.prettyOut("Create b0 Tmean")
            utils.Cmd( "fslmaths " +  os.path.join(dwi_path, "topup.dewarped_b0") + " -Tmean " + os.path.join(dwi_path, "topup.dewarped_b0.Tmean") + " -odt short").run()
            
        #    mask the average susceptibility distortion corrected b0 mean volume
        utils.Cmd("3dAutomask -dilate 2 -peels 0 -prefix " + mask_file + " " + os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz") ).run()
        
        # fsl might complain after above step because pixdim3 is 1.4999999 instead of 1.5
        # fslcpgeom topup.dewarped_b0.Tmean.nii.gz sub-83711_ses-1_dwi-edc_mask.nii.gz
        utils.Cmd("fslcpgeom " + os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz" ) + " " + mask_file).run()
        
        
    # snapshot mask
    o_file = os.path.join(subj.derivatives_path('dwi'), 'qc', 'mask')
    if not os.path.isfile(o_file + '.sag.png'):
        print('\n   +++ create mask qc images\n')
        u_lay = os.path.join(dwi_path, "topup.dewarped_b0.Tmean.nii.gz")
        o_lay = mask_file
        utils.Cmd("@chauffeur_afni -blowup 4 -ulay {} -olay {} -set_xhairs OFF -montx 3 -monty 1 -delta_slices 1 1 1 -set_subbricks 0 0 0 -ulay_range 0 5000 -save_ftype PNG -prefix {}".format(u_lay, o_lay, o_file)).run()





def get_highest_image_multiple(image_file_path):
    """
    Used to determine cnf file for FSL
    """
    if os.path.isfile(image_file_path):
        try:
            import nibabel as nib
            img = nib.load(image_file_path)
            
            highest_multiple = 1
            for mod in [4,2]:
                if ((img.header["dim"][1] % mod) + (img.header["dim"][2] % mod) + (img.header["dim"][3] % mod)) == 0:
                    highest_multiple = mod
                    return str(highest_multiple)
 
            return str(highest_multiple)
        
        except:
            raise            
    else:
        raise utils.CmdProcError("ERROR: get_highest_image_multiple() input not found")
        
        
        
        
        
        
def get_json_file_list(subj, indx=None, substr=""):
    """
    Return list of json sidecar files for subject visit
    If no files, returns empty list
    """
    
    import glob
    dwi_raw_path = subj.bids_modality_path(bidsType='dwi')
    
    # get dwi json sidecar
    #json_file_list = glob.glob(os.path.join(subj.bids_modality_path(bidsType='dwi'), "*.json"))
    json_file_list = sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.json'), key = lambda x: x.split('_')[-2])
    if not len(json_file_list):
        raise ValueError("   +++ Could not locate json sidecar at:\n    +++ "+ dwi_raw_path)
        
    if indx == None:
        if substr == "":
            #return glob.glob(os.path.join(subj.bids_modality_path(bidsType='dwi'), "*.json"))
            return sorted(glob.glob(dwi_raw_path + os.path.sep + '*run*_dwi.json'), key = lambda x: x.split('_')[-2])
        else:
            return glob.glob(os.path.join(subj.bids_modality_path(bidsType='dwi'), "*"+substr+"*.json"))
    elif ( isinstance(indx, int) and (len(json_file_list) > int(indx)) ):
        return json_file_list[indx]
    else:
        raise ValueError("   +++ Could not locate json sidecar at:\n    +++ "+ dwi_raw_path)
    
    
    
    
    
def get_json_param_value_dict(json_file_path, params=[]):
    """
    Return dictionary of requested keys with their corresponding values
    from the json file. If the key does not exist in the json file it will
    not be in the return dictionary
    """
    
    import json
    value_dict = {}
    with open(json_file_path) as fp:
        data = json.load(fp)
        for key in params:
            if key in data.keys():
                value_dict[key] = data[key]
    
    return value_dict
        
        
        


def clone_json_file_params(json_file_in_path, json_file_out_path, params=[]):
    """
    Copy a subset of json file parameters to a new json file
    """
    
    import json

    try:
        
        out_dict = {}
        with open(json_file_in_path, 'r') as f_in:
            data_in = json.load(f_in)
            for key in params:
                if key in data_in.keys():
                    out_dict[key] = data_in[key]
                else:
                    utils.warningMsg('key {} not found in json source file. Skipping.'.format(key))
                    
        with open(json_file_out_path, 'w') as f_out:
            json.dump(out_dict, f_out, indent = 6)
        
        print("   +++ Created:", json_file_out_path)

    except:
        raise 






def create_params_files(subj):
    """
    Use json sidecar input to create parameter files used by topup, applytopup, and eddy
    # n.b. Other parts of this pipeline currently assume the following:
    # dwi-params.txt contains only two lines, one for AP and the other for PA
    # dwi-index.txt has only "1"s and "2"s
    # this code was started so that the pipeline could be run in a container.
    # Don't use this unless you are putting dwi pipeline in a container and you change how this works
    
    """

    params_file = os.path.join(subj.derivatives_path(bidsType='dwi'), subj.bids_file_join('dwi-params.txt') )
    if not os.path.isfile(params_file):
        
        from mayerlab.preprocessing import bids

        # get json files in order
        dwi_json_list = get_json_file_list(subj)

        # write a line to params.txt for each json file
        with open(params_file, 'w') as f_out:
            for json_file in dwi_json_list:
                
                dwi_param = {'i':'0', 'j':'0', 'k':'0', 'TotalReadoutTime':'0'}
                pe = bids.SideCar.get_param(json_file, "PhaseEncodingDirection")
                param_key = pe[0]
                if len(pe) == 2:
                    param_val = '-1'
                else:
                    param_val = '1'
                    
                dwi_param[param_key] = param_val
                TotalReadoutTime = str(bids.SideCar.get_param(json_file, "TotalReadoutTime"))
                dwi_param["TotalReadoutTime"] = TotalReadoutTime
                f_out.write(" ".join(dwi_param.values()) + '\n')
                
        print("   +++ Created: " + params_file)
            


    index_file = os.path.join(subj.derivatives_path(bidsType='dwi'), subj.bids_file_join('dwi-index.txt') )
    if not os.path.isfile(index_file):
        
        from mayerlab.preprocessing import bids

        # get json files in order
        dwi_json_list = get_json_file_list(subj)
        
        # write a line to index.txt for each json file
        with open(index_file, 'w') as f_out:
            for file_indx, json_file in enumerate(dwi_json_list):
                bval_file = json_file.replace('json', 'bval')
                with open(bval_file, 'r') as bval_in:
                    v = len(bval_in.read().split())
                    
                indx_str = str(file_indx+1) + '\n'
                f_out.write( indx_str * v)
                
        print("   +++ Created: " + index_file)
        
        
        

        

    
class model:
    """
    models fit to DWI data
    """
    
    def generate_nonlinear_DT(subj, use_bvals='all'):
    
        """
        Standard AFNI nonlinear diffusion tensor fit.
    
        :param subj: utils.Subj
        :param study: utils.Study
        :param use_bvals: which gradients to use [some b-value,'all']
    
        """
    
        dwi_path = subj.derivatives_path(bidsType='dwi')
        
        in_bvals_path = os.path.join(dwi_path, subj.bids_file_join('dwi-concat.bvals-rounded'))
        if not os.path.isfile(in_bvals_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvals_path)
    
        in_dwi_path = os.path.join(dwi_path, subj.bids_file_join('dwi-edc.nii.gz'))
        if not os.path.isfile(in_dwi_path):
            raise utils.CmdProcError("Missing input:\n"+in_dwi_path)
    
        in_bvecs_path = os.path.join(dwi_path, subj.bids_file_join('dwi-edc.eddy_rotated_bvecs'))
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)
    
        in_mask_path = os.path.join(dwi_path, subj.bids_file_join('dwi-edc_mask.nii.gz'))
        if not os.path.isfile(in_mask_path):
            raise utils.CmdProcError("Missing input:\n"+in_mask_path)
            
            
        # dti output path
        out_dti_path = os.path.join(dwi_path ,"dti")
        utils.make_path(out_dti_path)                
        
        # just clean up name if using all b values.
        if use_bvals != 'all':
            out_dti_pre = os.path.join(out_dti_path, subj.bids_file_join("DTI_NL.bvecs-"+use_bvals+".nii.gz"))
            o_file_test = os.path.join(out_dti_path, subj.bids_file_join("DTI_NL.bvecs-"+use_bvals+"_FA.nii.gz"))
        else:
            out_dti_pre = os.path.join(out_dti_path, subj.bids_file_join("DTI_NL.nii.gz"))
            o_file_test = os.path.join(out_dti_path, subj.bids_file_join("DTI_NL_FA.nii.gz"))
        
        
        
        if not os.path.isfile(o_file_test):
            
            # by default use all gradients but if some b-value list is passed, use just those
            if use_bvals == 'all':
        
                brik_selector_str = ''
                #import shutil
                #shutil.copy(in_bvecs_path, os.path.join(out_dti_path, "3dDWItoDT.bvecs."+use_bvals))
                utils.Cmd('1dtranspose {} > {}'.format(in_bvecs_path, os.path.join(out_dti_path, "3dDWItoDT.bvecs."+use_bvals))).run()
        
            else:
        
                select_bvals = [0]
                select_bvals.append(int(use_bvals))
                brik_selector_ind = []
        
                brik_indx = 0
                with open(in_bvals_path) as f:
                    for line in f:
                        for bval in line.strip().split(' '):
        
                            bval_i = int(round(float(bval+".0")/100)*100)
                            if bval_i in select_bvals:
                                brik_selector_ind.append(brik_indx)
        
                            brik_indx += 1
        
        
                brik_selector_str = "'["+ ','.join(map(str, brik_selector_ind)) + "]'"
        
        
        
                new_bvecs_file = open(os.path.join(out_dti_path, "3dDWItoDT.bvecs."+use_bvals), "w")
                with open(in_bvecs_path) as f:
                    brik_indx = 0
                    for line in f:
                        if brik_indx in brik_selector_ind:
                            new_bvecs_file.write(line)
                        brik_indx += 1
        
        
                new_bvecs_file.close()
        
            

        
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Nonlinear Calculation of Tensor")
            utils.prettyOut("Using b="+use_bvals+" gradients for model")
    
            # create 3dDWItoDT version of bvecs by dropping initial b0 gradient
            utils.Cmd("tail -n +2 "+ os.path.join(out_dti_path, "3dDWItoDT.bvecs."+use_bvals) +" > " + os.path.join(out_dti_path, "3dDWItoDT.bvecs")).run()

            utils.Cmd("3dDWItoDT -prefix " + os.path.join(out_dti_path,out_dti_pre) +
                                " -mask " + in_mask_path + " " +
                                " -nonlinear -max_iter 10 -max_iter_rw 5 -reweight -eigs -sep_dsets -mean_b0 " +
                                " " + os.path.join(out_dti_path, "3dDWItoDT.bvecs") + " " +
                                in_dwi_path + brik_selector_str).run()
                
            utils.purge_paths( os.path.join(out_dti_path, "3dDWItoDT.bvecs*"), os.path.join(out_dti_path, "*DTI_NL*_DT.nii.gz") )
        




    
    
    
        
class qc:
    """
    functions that capture and output QC images for DWI
    """
        
        
    def snapshot_b0(subj, b0_filepath):
        
        # snapshot topup input pe 1
        o_file = os.path.join(subj.derivatives_path('dwi'), 'qc', 'topup_b0_input_0')
        if not os.path.isfile( o_file + '.sag.png'):
            print('\n   +++ snapshot topup input pe 1\n')
            utils.Cmd("@chauffeur_afni -blowup 4 -ulay {} -olay_off -set_xhairs OFF -montx 3 -monty 1 -delta_slices 1 1 1 -set_subbricks 0 0 0 -ulay_range 0 5000 -save_ftype PNG -prefix {}".format(b0_filepath, o_file)).run()
        
        # snapshot topup input pe 2
        o_file = os.path.join(subj.derivatives_path('dwi'), 'qc', 'topup_b0_input_1') 
        if not os.path.isfile( o_file + '.sag.png'):
            print('\n   +++ snapshot topup input pe 2\n')
            utils.Cmd("@chauffeur_afni -blowup 4 -ulay {} -olay_off -set_xhairs OFF -montx 3 -monty 1 -delta_slices 1 1 1 -set_subbricks 1 1 1 -ulay_range 0 5000 -save_ftype PNG -prefix {}".format(b0_filepath, o_file)).run()
        
        # snapshot topup output
        o_file = os.path.join(subj.derivatives_path('dwi'), 'qc', 'topup_b0_output')
        if not os.path.isfile( o_file + '.sag.png'):
            print('\n   +++ snapshot topup output\n')
            utils.Cmd("@chauffeur_afni -blowup 4 -ulay {} -olay_off -set_xhairs OFF -montx 3 -monty 1 -delta_slices 1 1 1 -set_subbricks 0 0 0 -ulay_range 0 5000 -save_ftype PNG -prefix {}".format(os.path.join(subj.derivatives_path('dwi'), 'topup.dewarped_b0.Tmean.nii.gz'), o_file)).run()


        
        
        
        
    def snapshot_spatial_normalization(subj, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
      
        # if not there, make it
        utils.make_path(os.path.join(subj.derivatives_path(bidsType='dwi'), 'qc'))
        
        if not os.path.isfile(os.path.join(subj.derivatives_path(bidsType='dwi'), 'qc', subj.bids_file_join('dwi', space+'.jpg'))):
            
            if space == 'MNI':
                apply_spatial_normalization(subj, 
                                                os.path.join(subj.derivatives_path('dwi'), 'topup.dewarped_b0.Tmean.nii.gz'),
                                                os.path.join(subj.derivatives_path('dwi'), 'qc', 'VERIFY.nii.gz'),
                                                space='MNI', template_key='mni_icbm152_nlin_asym_2009c')
                
                ulay_file = utils.Env.template_lib[template_key]['head']
                
            else:
                apply_spatial_normalization(subj, 
                                                os.path.join(subj.derivatives_path('dwi'), 'topup.dewarped_b0.Tmean.nii.gz'), 
                                                os.path.join(subj.derivatives_path('dwi'), 'qc', 'VERIFY.nii.gz'), 
                                                space='TLRC' )                

                ulay_file = utils.Env.template_lib['TT_N27']['head']
                
                
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : QC "+space+" Snapshots" )
            olay_file = os.path.join(subj.derivatives_path('dwi'), 'qc', 'VERIFY.nii.gz')
            utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, olay_file, os.path.join(subj.derivatives_path(bidsType='dwi'), 'qc', subj.bids_file_join('dwi', space)))).run()
            
            utils.purge_paths(os.path.join(subj.derivatives_path(bidsType='dwi'), 'qc', 'VERIFY*'))
    
            
    
    
        
        
    def aggregate_qc_pdf(subj, space='MNI'):
        """
        Create the Preprocessing QC PDF with a fixed set of images
        """    
        import pylatex
        from pylatex.utils import NoEscape        
        
        
        dwi_path = subj.derivatives_path(bidsType='dwi')
        
        o_file = os.path.join(dwi_path, 'qc', 'qc-preprocessing')
        if not os.path.isfile(o_file + '.pdf'):
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : QC PDF" )
            
            try:
                
                geometry_options = {"tmargin": "1cm", "lmargin": "3cm"}
                doc = pylatex.Document(page_numbers=False, geometry_options=geometry_options)
    
                doc.preamble.append(pylatex.Command('title', subj.study.label + ' DWI Preprocessing QC'))
                doc.preamble.append(pylatex.Command('author', subj.ursi.full + " : Visit " + str(subj.visit)))
                doc.preamble.append(pylatex.Command('date', NoEscape(r'\today')))
                doc.append(NoEscape(r'\maketitle'))
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('Spatial Normalization', numbering=False)):
                        doc.append('Colored edge overlay should align with tissue boundaries')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image( os.path.join(dwi_path, 'qc', subj.bids_file_join('dwi_'+space+'.jpg')), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('Distortion Corrected B0 registration over '+space+' template')
                        
                        
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('Distortion Correction', numbering=False)):
                        doc.append('Input images should contain minimal artifacts')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(dwi_path, 'qc', 'topup_b0_input_0.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('Uncorrected B0 AP Input - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(dwi_path, 'qc', 'topup_b0_input_1.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('Uncorrected B0 PA Input - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic3:
                            pic3.add_image(os.path.join(dwi_path, 'qc', 'topup_b0_output.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic3.add_caption('Corrected B0 Output - Sagittal View')   
                
                doc.append(NoEscape(r'\newpage'))
                
                with doc.create(pylatex.Center()):
                    with doc.create(pylatex.Section('DWI Mask', numbering=False)):
                        doc.append('Mask should cover entire brain')
                        with doc.create(pylatex.Figure(position='h')) as pic1:
                            pic1.add_image(os.path.join(dwi_path, 'qc', 'mask.sag.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic1.add_caption('DWI Mask - Sagittal View')
                        with doc.create(pylatex.Figure(position='h')) as pic2:
                            pic2.add_image(os.path.join(dwi_path, 'qc', 'mask.cor.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic2.add_caption('DWI Mask - Coronal View')
                        with doc.create(pylatex.Figure(position='h')) as pic3:
                            pic3.add_image(os.path.join(dwi_path, 'qc', 'mask.axi.png'), width=NoEscape('1.0\\linewidth'), placement=None )
                            pic3.add_caption('DWI Mask - Axial View')  
                        
                        
                doc.generate_pdf(o_file, clean_tex=True, compiler='pdflatex')


            except:
                
                # clean up attempt files on error
                utils.purge_paths(os.path.join(dwi_path, 'qc', 'qc-preprocessing.*'))
                raise

 
    
 
    
    
    def get_edc_image_outlier_summary(subj,study):
        """
        Return outlier summary from eddy_outlier_map, and qc
        """
        # dwi_path = subj.path_append(study, proc_dir)
        # base_name   = study.label + "." + subj.ursi.short + "." + "dwi.concat.edc"
        
        # import numpy as np
        
        # # load outlier matrix
        # try:
            
        #     i_file = os.path.join(dwi_path, base_name + '.eddy_outlier_map')
        #     if os.path.isfile(i_file):
        #         outlier_map = np.loadtxt( os.path.join(dwi_path, base_name + '.eddy_outlier_map'), skiprows=1, dtype='int')
        #         total_num_slices = outlier_map.shape[0] * outlier_map.shape[1]
        
        #         sum = np.sum(outlier_map, axis=1)
        #         pct = np.round((sum/outlier_map.shape[1])*100)
        #         r = np.array([list(range(outlier_map.shape[0])), sum, pct.astype(dtype='int')])
                
        #         outlier_slice_sum =np.sum(r[1])
        #         outlier_image_count1 = np.count_nonzero(r[1])
        #         outlier_image_count3 = np.count_nonzero(r[1]>3)
        #         outlier_slice_pct = (outlier_slice_sum/total_num_slices)*100
        #         outlier_image_count1_pct = (outlier_image_count1/outlier_map.shape[0])*100
        #         outlier_image_count3_pct = (outlier_image_count3/outlier_map.shape[0])*100
                
        #         return outlier_slice_sum, outlier_slice_pct, outlier_image_count1, outlier_image_count1_pct, outlier_image_count3, outlier_image_count3_pct, total_num_slices
                
        #     else:
        #         print("{} couldn't be located.".format(i_file))
                
        # except:
        #     raise
            
            
    
    def view_edc_image_outlier_summary(subj,study):
        """
        Generate outlier summary from eddy_outlier_map
        """
        # # load outlier matrix
        # try:
        #     utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Gradient QC")
            
        #     ( outlier_slice_sum, outlier_slice_pct, outlier_image_count1, outlier_image_count1_pct, outlier_image_count3, outlier_image_count3_pct, total_num_slices ) = get_edc_image_outlier_summary(subj,study)
        #     print( "{} outlier artifacts found in {:0.2f}% of the slices.".format(outlier_slice_sum, outlier_slice_pct ) )
        #     print( "{} of the images({:0.2f}%) have 1 or more slice outliers. [Single Zipper].".format( outlier_image_count1, outlier_image_count1_pct ) )
        #     print( "{} of the images({:0.2f}%) have 3 or more slice outliers. [Multiple Zippers].".format( outlier_image_count3, outlier_image_count3_pct) )
                
        # except:
        #     raise
        
        
        
        
    def run_eddy_squad(study, visit="", subjects=[], output_dir=""):
        """
        Generate study-wise qc from eddy_quad output
        
        """
    
        # input_list = []
        
        # for subject_id in subjects:
        #     subj = utils.Subj(subject_id, visit=visit)
        #     qc_path = subj.path_append(study, proc_dir, 'qc')
        #     if os.path.isdir( qc_path ):
        #         input_list.append(qc_path)
                
        # utils.Cmd("eddy_squad " + " ".join(input_list) + " --update --output-dir="+output_dir).run()
        
        
        
        
    def view_eddy_quad_qc_stats(subj, study, raw=False):
        """
        Returns summary stats from quad report
        """
        
        # try:
        #     qc_data = get_eddy_quad_qc_json(subj, study)
            
        #     utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : EDDY QC Report Stats")
        #     print("Total Outliers (%): {:0.2f}".format(qc_data["qc_outliers_tot"]))
        #     print("Volume-to-Volume Motion:")
        #     print("   Average Abs motion (mm): {:0.2f}".format(qc_data["qc_mot_abs"]))
        #     print("   Average Rel motion (mm): {:0.2f}".format(qc_data["qc_mot_rel"]))
            
        # except:
        #     raise




    def run_eddy_quad(subj):
        """
        Generate eddy qc output
        """
        
        import shutil
        
        dwi_path = subj.derivatives_path(bidsType='dwi')
        
        edc_basename_path = os.path.join(dwi_path, subj.bids_file_join("dwi-edc"))
        mask_path = os.path.join(dwi_path, subj.bids_file_join('dwi-edc_mask.nii.gz'))
        acqp = os.path.join(subj.study.paramsPath, 'dwi_params.txt')
        index = os.path.join(subj.study.paramsPath, 'dwi_index.txt')
        slspec = os.path.join(dwi_path, 'dwi.slspec' )
        # get dwi json sidecar for slicetiming
        json_file = get_json_file_list(subj, indx=0)
                    
        
        try:
            if not os.path.isfile( os.path.join(dwi_path, 'qc','qc.json' )):
        
                # eddy_quad won't proceed if directory already exits. i will not assure that
                # some other qc operation hasn't already created it so we will work around
                # eddy_quad by sending output to temp directory and then copying that to qc
                
                
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Running QC")
                
                if not os.path.isfile(slspec):
                    generate_fsl_eddy_slspec(json_file, slspec)
                
                # go into dwi directory because of the way eddy_quad uses the full path in display
                os.chdir(dwi_path)
                utils.make_path(os.path.join(dwi_path, 'qc'))
                
                utils.Cmd("eddy_quad " + edc_basename_path
                          + " -v --eddyIdx " + index
                          + " --eddyParams " + acqp
                          + " --mask " + mask_path
                          + " --bvals " + os.path.join(dwi_path, subj.bids_file_join("dwi-concat.bvals"))
                          + " --bvecs " + os.path.join(dwi_path, subj.bids_file_join("dwi-edc.eddy_rotated_bvecs"))
                          + " --slspec " + slspec
                          + " --output-dir eddy_quad_tmp " ).run()
                
                
                # rename pdf
                os.rename( os.path.join(dwi_path, 'eddy_quad_tmp', 'qc.pdf'), os.path.join(dwi_path, 'eddy_quad_tmp', 'qc-eddy.pdf') )
                
                # no move contents to qc
                for f in os.listdir(os.path.join(dwi_path, 'eddy_quad_tmp')):
                    shutil.move(os.path.join(dwi_path, 'eddy_quad_tmp', f), os.path.join(dwi_path, 'qc/'))
                    

                
        except:
            raise
            
        finally:
            utils.purge_paths(os.path.join(dwi_path, 'eddy_quad_tmp'))
    
    
    
    
    
    def get_eddy_quad_json(subj, study):
        """
        Returns stats from quad report in json format
        """
        
        try:

            dwi_path = subj.derivatives_path(bidsType='dwi')
            qc_file = os.path.join(dwi_path, 'qc','qc.json')
            if not os.path.isfile(qc_file):
                qc.run_eddy_quad(subj,study)
        
            # get stats
            import json
    
            with open(qc_file,'r') as f:
                qc_data = json.load(f)
                return qc_data
    
        except:
            raise






class denoise:
    """
    functions for denoising DWI
    """
    
    def MPPCA(subj):
        """
        dMRI noise level estimation and denoising using Marchenko-Pastur PCA
        
        DWI data denoising and noise map estimation by exploiting data redundancy in the PCA domain using the prior knowledge that the eigenspectrum of random covariance matrices is described by the universal Marchenko-Pastur (MP) distribution. Fitting the MP distribution to the spectrum of patch-wise signal matrices hence provides an estimator of the noise level ‘sigma’, as was first shown in Veraart et al. (2016) and later improved in Cordero-Grande et al. (2019). This noise level estimate then determines the optimal cut-off for PCA denoising.

        Important note: image denoising must be performed as the first step of the image processing pipeline. The routine will fail if interpolation or smoothing has been applied to the data prior to denoising.

        Note that this function does not correct for non-Gaussian noise biases present in magnitude-reconstructed MRI images. If available, including the MRI phase data can reduce such non-Gaussian biases, and the command now supports complex input data.
        """

        
        dwi_derv_path = subj.derivatives_path(bidsType='dwi')
        i_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi-concat.nii.gz"))
        o_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi-concat_den-mppca.nii.gz"))
        noise_map_file = os.path.join(dwi_derv_path, subj.bids_file_join("dwi_den-noise-map.nii.gz"))
        
        if not os.path.isfile(noise_map_file):
            try:
                utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : DWI : Denoise")
                utils.Cmd("dwidenoise -info -estimator Exp2 -nthreads 8 -noise " + noise_map_file + " "  + i_file + " " + o_file).run()
                # float32 output makes no sense, turn it back into int32
                utils.Cmd("fslmaths " + o_file + " " + o_file + " -odt int").run()

                # since first step in the processing we'll just replace concat file with denoised file
                utils.purge_paths(i_file)
                os.rename(o_file, i_file)
                
            except:
                raise        

        
        
        
    def nlmeans(subj,study):
    
        """
        Denoises DWI images using Non-Local Means (NLMEANS)
    
        :param subj: utils.Subj
        :param study: utils.Study
    
        """
    
        import nibabel as nib
        from dipy.denoise.nlmeans import nlmeans
        from dipy.denoise import noise_estimate
    
    
        dwi_path = subj.derivatives_path(bidsType='dwi')
        in_file = os.path.join(dwi_path, subj.bids_file_join("dwi-concat.nii.gz"))
        o_file = os.path.join(dwi_path, subj.bids_file_join("dwi-concat-denoised.rician.est_sigma.nii.gz"))
        
        if not os.path.isfile(o_file):
            #    load image
            img = nib.load(in_file)
            aff = img.get_affine()
            #    pull data from image object
            data = img.get_data()
        
            #    noise estimate per image ( Standard deviation estimation from local patches http://nipy.org/dipy/reference/dipy.denoise.html )
            est_sigma = noise_estimate.estimate_sigma(data, disable_background_masking=False)
            #    calculate denoised dwi
            den = nlmeans(data, sigma=est_sigma, mask=None, patch_radius=1, block_radius=5, rician=True)
            #    write out denoised volume
            nib.save(nib.Nifti1Image(den, aff), o_file)
    
