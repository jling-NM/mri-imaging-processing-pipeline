# -*- coding: utf-8 -*-

"""@package asl
Preprocessing module for arterial spin label imaging

Created on Wed Dec  2 13:28:49 2015

@author: jling

"""







def generate_3d_pcasl_cbf(i_path, i_file_pcasl, i_file_m0, o_file_cbf, i_file_mask=None, o_path=None, pld=1800, m0_scale=10):

    """    
    
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
        tau      (label duration of the labeling RF pulse train)     [xxx]

    
        #        Tisdall, Dylan <mtisdall@pennmedicine.upenn.edu>
        #        Sep 20, 2018, 9:21 PM (4 days ago) 	
        #        to John, Andy, Arvind, me
        #        Hi Andy,
        #        
        #        Just chiming in to answer those two questions about scaling:
        #        > a) What if any scaling is present in the perfusion image calculated on the scanner? It does not seem to be purely the mean of the difference images? 
        #        I'm not 100% certain of all the steps that go into computing the MeanPerf image, but my understanding is that it is at least scaled to cover the range 0-4096 and truncated below at 0. As such, I wouldn't expect it to match just the mean label-control difference exactly. As John mentioned, I think the MeanPerf series is best used as an on-scanner diagnostic that you're getting reasonable data.
        #        > b) When we change dicom factor (currently we -1500) does it equally effect the ASL and M0  images?  There is another M0 factor set to 10 - does that matter in the CBF calculations?  We are ignoring it.       
        #        Yes, the DICOM scale factor affects both the ASL and M0 images. The M0 factor being set to 10 means that the M0 data is scaled like the ASL data, and then additionally divided by 10 before being saved as a DICOM. This is to ensure its dynamic range is roughly the same as the ASL images, and so will be well-scaled for the DICOM images. When you're doing your analysis, you'll need to multiply the M0 value at each voxel by 10 to get the true value of M0. If you don't do this, your estimates of the \DeltaM / M0 will be too big by a factor of 10.

        #        Cheers,
        #        Dylan


    """

    import os
    import numpy as np
    import numpy.ma as ma
    from mayerlab import utils
    from termcolor import cprint    #    for colored printing


    try:

        ##################################################
        # OPEN IMAGES
        ##################################################

        print("Loading data...")

        # load asl mask volume first
        if i_file_mask is not None:
            mask_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_mask))

        # load pcasl volume
        asl_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_pcasl), dtype='int16')
        # unpack dimensions to make code more readable
        (asl_xdim,asl_ydim,asl_zdim,asl_tdim) = asl_AfniImg.data.shape

        
        if i_file_m0 is not None:
            # load perf volume
            m0_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_m0), dtype='int16')
            
            # when i create m0 vol i'm do so as a masked array so that when we use it in division.
            if i_file_mask is not None:
                m0_vol = m0_scale * ma.masked_where(mask_AfniImg.data==0, m0_AfniImg.data, copy=True)
            else:
                m0_vol = m0_scale * m0_AfniImg.data



        ##################################################
        # DEFINE VALUES FOR CBF CALCULATION
        ##################################################

        # blood/tissue water partition coefficient
        # 0.9 for pCASL Alsop 2014
        # [mL/g]
        blood_tissue_h2o_part_coef = 0.9

        # label_duration (aka Labeltime) =  Number of RF Blocks*rf_block_duration;
        #label_duration = (num_rf_blocks*(rf_block_duration*1000))
        label_duration = 1800

        # r1a=1/BloodT1; %0.601sec for 3T
        # the longitudinal relaxation time of blood
        # 1650 for pCASL Alsop 2014
        # [msec]
        T1_blood = 1650
        
        # labeling efficiency, 0.85 for pCASL Alsop 2014
        #labeling_efficiency = 0.85
        
        #        Detre, John <detre@pennmedicine.upenn.edu>
        #        AttachmentsSep 20, 2018, 3:14 PM (4 days ago)
        #        
        #        Hi Andy,
        #        
        #        With regard to GM:WM resolution, there is considerable Z-blurring in 3D imaging that will lower this ratio.  Kz acceleration helps with that, but when we want to isolate GM and WM we use a 4-shot protocol with 2.5 mm isotropic resolution.
        #        Also, background suppression is known to reduce labeling efficiency.  We use 0.72 instead of 0.85.
        #        John
        labeling_efficiency = 0.72


        print("\n")
        cprint("\n=======================================================", 'red', 'on_grey', attrs=['bold'] )
        cprint("  CBF Parameters", 'red', 'on_grey', attrs=['bold'] )
        cprint('  + blood_tissue_h2o_part_coef : {0:.2f}'.format(blood_tissue_h2o_part_coef), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + pld : {0:d}'.format(pld), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + m0_scale : {0:d}'.format(m0_scale), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + label_duration : {0:.2f}'.format(label_duration), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + T1_blood : {0:.2f}'.format(T1_blood), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + labeling_efficiency : {0:.2f}'.format(labeling_efficiency), 'red', 'on_grey', attrs=['bold'] )
        if i_file_m0 is not None:
            cprint('  + using PD image for m0 normalization', 'red', 'on_grey', attrs=['bold'] )
        else:
            cprint('  + using control image normalization', 'red', 'on_grey', attrs=['bold'] )
        cprint("=======================================================", 'red', 'on_grey', attrs=['bold'] )

        

        # perfusion image is mean of all labeled images subtracted from all control or unlabeled images (control - label)
        mean_diff = np.mean(asl_AfniImg.data[:,:,:,1::2], axis=3) - np.mean(asl_AfniImg.data[:,:,:,::2], axis=3)
        
        
        # calculate mean cbf using alsop formula
        #          6000 * λ * (Control-Tagged) * e^(PLD/T1_blood)
        #    CBF = ---------------------------------------------------  [ml/100g/min]
        #          2 * alpha * T1_blood * PD * (1 - e^-(tau/T1_blood))        

        # The factor of 6000 converts the units from ml/g/s to ml/(100g)/min, which is customary in the physiological literature.
        # The 1000 converts from milleseconds to seconds. The input parameters are given in milleseconds in Alsop 2014.
        # To help in reporting i've left the inputs in milleseconds and the two rescaling factors separate in the formula. i.e. inputs in milleseconds, output in minutes

        # divide by zero errors could be a good diagnostic but i have found it to be driven by non-parenchyma while steps below clean up the debris
        # therefore, turn off these errors for the calculation
        with np.errstate(divide='ignore', invalid='ignore'):
            mean_cbf_vol = ( 6000 * 1000 * blood_tissue_h2o_part_coef * mean_diff * np.exp(pld/T1_blood) ) / ( 2 * labeling_efficiency * T1_blood * m0_vol * (1-(np.exp(-(label_duration/T1_blood)))) )
            
        # convert INF to 0 to forestall conversion to +/-1.7976931348623157e+308 by np.nan_to_num().
        mean_cbf_vol[np.isinf(mean_cbf_vol)] = 0.0
        
        # remove nan's that might occur above from divide-by-zero
        # this will also convert INF to +/-1.7976931348623157e+308 which i don't find useful as it won't store in float32
        # and is a useless arfifact of divide by zero
        mean_cbf_vol = np.nan_to_num(mean_cbf_vol)
        
        # finally, write out mean cbf volume
        print("Writing Mean CBF Volume...")

        if o_path is None:
            o_path = i_path

        # create output image from input image
        if i_file_m0 is not None:
            cbf_img = utils.AfniImage(np.array(mean_cbf_vol,dtype=np.float32), m0_AfniImg.header)
        else:
            cbf_img = utils.AfniImage(np.array(mean_cbf_vol,dtype=np.float32), asl_AfniImg.header)
        
        # because save function automatically adds space let's remove it if it shows up here.
        o_file_cbf = o_file_cbf.replace('+orig', '')
        
        # save the image to file
        cbf_img.save(os.path.join(o_path, o_file_cbf), gz=True)


    except:

        raise




def generate_pcasl_cbf(i_path, i_file_pcasl, i_file_m0, o_file_cbf, i_file_mask=None, o_path=None, min_tr_ms=4140, pld=1800, T1_blood=1650):

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

     





    i_path
    Reguired
    From where is data coming

    i_file_pcasl
    Reguired
    The pCASL input
    
    i_file_m0
    Reguired
    m0 data to use for normalization
    
    o_prfx
    Reguired
    Attach this to beginning of CBF output file name

    i_file_mask
    Optional
    Mask the calculations and CBF output

    o_path
    Required
    Where is output to go

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
    from mayerlab import utils
    from termcolor import cprint    #    for colored printing


    try:

        ##################################################
        # OPEN IMAGES
        ##################################################

        print("Loading data...")

        # load asl mask volume first
        if i_file_mask is not None:
            mask_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_mask))

        # load pcasl volume
        asl_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_pcasl), dtype='int16')
        # unpack dimensions to make code more readable
        (asl_xdim,asl_ydim,asl_zdim,asl_tdim) = asl_AfniImg.data.shape

        
        if i_file_m0 is not None:
            # load perf volume
            m0_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_m0), dtype='int16')
            
            # when i create m0 vol i'm do so as a masked array so that when we use it in division.
            if i_file_mask is not None:
                m0_vol = ma.masked_where(mask_AfniImg.data==0, m0_AfniImg.data, copy=True)
            else:
                m0_vol = m0_AfniImg.data



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
        if i_file_m0 is not None:
            cprint('  + using PD image for m0 normalization', 'red', 'on_grey', attrs=['bold'] )
        else:
            cprint('  + using control image normalization', 'red', 'on_grey', attrs=['bold'] )
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
        mean_diff = np.mean(asl_AfniImg.data[:,:,:,1::2], axis=3) - np.mean(asl_AfniImg.data[:,:,:,::2], axis=3)


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
        mean_cbf_vol = np.nan_to_num(mean_cbf_vol)
        
        masked_data = mean_cbf_vol.data * mask_AfniImg.data
        # remove negative CBF which makes no sense and implausibly high CBF. oxford_asl also zeros out negative values
        masked_filtered_data = np.where( (masked_data < 0.0) | (masked_data > 200.0), 0, masked_data)
                
                
        # finally, write out mean cbf volume
        print("Writing Mean CBF Volume...")

        if o_path is None:
            o_path = i_path

        # create output image from input image
        if i_file_m0 is not None:
            cbf_img = utils.AfniImage(np.array(masked_filtered_data, dtype=np.float32), m0_AfniImg.header)
        else:
            cbf_img = utils.AfniImage(np.array(masked_filtered_data, dtype=np.float32), asl_AfniImg.header)
        
        # because save function automatically adds space let's remove it if it shows up here.
        o_file_cbf = o_file_cbf.replace('+orig', '')
        
        # save the image to file
        cbf_img.save(os.path.join(o_path, o_file_cbf), gz=True)

    except:

        raise





def qa_cbf(subj, study):
    
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

    from mayerlab import utils
    import os
    import seaborn as sns
    import matplotlib.pyplot as plt

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : CBF distributions" )

    os.chdir("/tmp")
    
    # move cbf native to subject's T1 space
    i_file_cbf = os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.mean.cbf.fwhm6+orig.HEAD")
    o_file_cbf = "/tmp/VERIFY.cbfT1+orig.HEAD"
    if not os.path.isfile(i_file_cbf):
        raise utils.CmdProcError("Missing: "+i_file_cbf)
        
    if not os.path.isfile(o_file_cbf):
        utils.Cmd("3dAllineate " +
                " -master " + subj.get_path(study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.sksp+orig " +
                " -1Dmatrix_apply " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + ".pCASL_epi2T1.aff12.1D" +
                " -input " + i_file_cbf +
                " -floatize" +
                " -interp trilinear" +
                " -final wsinc5"+
                " -prefix " + o_file_cbf ).run()

    
    # load wm mask
    i_file_wm = subj.get_path(study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm.mask.60+orig.HEAD"
    if not os.path.isfile(i_file_wm):
        raise utils.CmdProcError("Missing: "+i_file_wm)
    img_wm = utils.AfniImage.load(i_file_wm)
    wm_mask = img_wm.data
    
    # load gm mask
    i_file_gm = subj.get_path(study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.spm.gm.mask.60+orig.HEAD"
    if not os.path.isfile(i_file_wm):
        raise utils.CmdProcError("Missing: "+i_file_gm)
    img_gm = utils.AfniImage.load(i_file_gm)
    gm_mask = img_gm.data
    
    # load T1 space CBF
    img_cbf = utils.AfniImage.load(o_file_cbf)
    cbf_data = img_cbf.data
    
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
    # display plots
    plt.show()
    plt.close()

    utils.purge_paths("/tmp/VERIFY.cbfT1+orig*")
    
    #plt.title("CBF masked distributions")
    #plt.xlabel("CBF")
    # plot masked distributions
    i_file_fd = os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.3dreg.log.raw.deriv.fd.rowsum")

    with open(i_file_fd) as f:
        fd_data = f.read()
    
    fd_data = list(map(float, fd_data.replace(' ', '').split('\n')[:-1]))
    plt.plot(fd_data)
    plt.title("PCASL Framewise Displacement")
    plt.xlim([0,90])
    plt.xlabel("Image")
    plt.ylabel("FD")
    sns.despine(top=True, right=True, left=True, bottom=False, offset=None, trim=False)
#    sns.kdeplot( cbf_data[gm_mask.nonzero()].flatten(), shade=False, label="CBF in GM", linewidth=1.5 );
#    # remove top/right axis lines
#    sns.despine(top=True, right=True, left=True, bottom=False, offset=None, trim=False)
#    plt.xlim(-50.0, 100.0)
#    # display plots
    plt.show()
#    plt.close()
    
    
    
    


def qa_cbf_mask(subj, study, mask_path):
    
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
    i_file_cbf = os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.mean.cbf.fwhm6+tlrc.HEAD")
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

    
    
    
def generate_pcasl_cbf_wang(i_path, i_file_pcasl, o_prfx, i_file_mask=None, o_path=None, i_file_m0=None, min_tr_ms=4140, tag_delay_ms=1800):

    """

    i_path
    Reguired
    From where is data coming

    i_file_pcasl
    Reguired
    The pCASL input

    o_prfx
    Reguired
    Attach this to beginning of CBF output
    
    i_file_mask
    Optional
    Mask the calculations and CBF output

    o_path
    Required
    Where is data going

    i_file_m0
    Optional
    m0 data to use for normalization. If not provided, the pair control image used instead

    min_tr_ms
    Required
    See inline comments in code

    tag_delay_ms
    Required
    See inline comments in code




     CBF quantification:
     For casl, the CBF value is calculated by:
       CBF (ml/100g/min) = 60*100*deltaM*lambda*R/(2*alpha*M0*(exp(-w*R)-exp(-(tao+w)*R))
       where deltaM - ASL perfusion difference, lambda - blood/tissue water partition coefficient, R - longitudinal relaxation rate of blood,
           alpha - tagging efficiency,
           M0 - equilibrium magnetization of brain. M0 should be acquired when the spins are fully relaxed (with long TR and short TE)
           w - post-labeling delay, w is generally assumed to be longer than the transit time from the tagging place to the arterial vascular compartment
               to allow the delivery of all the tagged blood into the imaging slices. See Alsop and Detre, J Cereb Blood Flow Metab 1996;16:1236每1249.
           tao - duration of the labeling RF pulse train.
       regular values for some parameters are lambda=0.9g/ml,
       for 3T, alp=0.68, R=1/1664ms  was R=0.67sec-1 in Wang 03  but should be smaller than that according to the literature.
       for 1.5T, alp=0.71, R=0.83sec-1.
       Please refer to: Wang J, Alsop DC, et al. (2003) Arterial transit time imaging with flow encoding arterial spin tagging (FEAST).
       Magn Reson Med. 50(3), page600, formula [1]
         might be changed to
       CBF (ml/100g/min) = 60*100*deltaM*lambda*R1tissue/(2*alpha*M0*(exp(-w*R1b)(1-exp(-(tao)*R1tissue)))  according to Buxton's model.
     - FROM ASLTOOLBOX


     Compare with Alsop et al. 2014:

              6000 * lambda * (Control-Tagged) * e^(PLD/T1_blood)
       CBF = ---------------------------------------------------
              2 * alpha * T1_blood * PD * (1 - e^-(tau/T1_blood))

        Where:
            lambda   = brain/blood partition coef (mL/g)            [0.9 mL/g]
            PLD      = post-label delay
            T1_blood = longitudinal relaxation time of blood (secs) [1650 ms]
            alpha    = labelling efficiency                         [0.85]
            PD       = spin-free scaling image
            tau      = label duration

    """

    import os
    import numpy as np
    import numpy.ma as ma
    from mayerlab import utils
    from termcolor import cprint    #    for colored printing


    try:

        ##################################################
        # OPEN IMAGES
        ##################################################

        print("Loading data...")

        # load asl mask volume first
        if i_file_mask is not None:
            mask_AfniImg = utils.AfniImage.load(os.path.join(i_path,i_file_mask))

        # load pcasl volume
        asl_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_pcasl), dtype='int16')
        # unpack dimensions to make code more readable
        (asl_xdim,asl_ydim,asl_zdim,asl_tdim) = asl_AfniImg.data.shape
        # get a count of the number of imaging pairs i.e. label-control
        img_pairs = asl_tdim/2

        if i_file_m0 is not None:
            # load perf volume
            m0_AfniImg = utils.AfniImage.load(os.path.join(i_path, i_file_m0), dtype='int16')
            
            # when i create m0 vol i'm do so as a masked array so that when we use it in division later it does not throw errors
            if i_file_mask is not None:
                m0_vol = ma.masked_where(mask_AfniImg.data==0, m0_AfniImg.data, copy=True)
            else:
                # JL : 20180531 : removed masking. Experimental
                m0_vol = m0_AfniImg.data



        ##################################################
        # DEFINE VALUES FOR CBF CALCULATION
        ##################################################

        # blood/tissue water partition coefficient
        # [mL/g]
        blood_tissue_h2o_part_coef = 0.9

        # keep clear what these are in calculations
        num_rf_blocks = asl_tdim
        num_slices = asl_zdim

        # tagging_time (aka Labeltime) =  Number of RF Blocks*0.0185;
        # RORDEN: "The CFN pCASL RF block duration is ALWAYS = 0.0185s (20 RF pulses with gaps) - 18500us"
        rf_block_duration = 0.0185
        # Number of RF Blocks = "Num RF Blocks" from sequence protocol  e.g. 82*0.0185*1000
        # [msec]
        tagging_time = int(num_rf_blocks*rf_block_duration*1000)


        # Post tagging delay/Post-labeling delay; get from sequence protocol
        # [msec]
        post_tag_delay_time = tag_delay_ms

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
        slice_time = (min_tr_ms-tagging_time-post_tag_delay_time)/float(num_slices)


        # r1a=1/BloodT1; %0.601sec for 3T
        # the longitudinal relaxation time of blood
        # [msec]
        r1a = 0.00060096

        # labeling efficiency, 0.85 for pCASL from Dai,W,Alsop - Magn Reson Med 2008
        tagging_efficiency = 0.85


        print("\n")
        cprint("\n=======================================================", 'red', 'on_grey', attrs=['bold'] )
        cprint("  CBF Parameters", 'red', 'on_grey', attrs=['bold'] )
        cprint('  + blood_tissue_h2o_part_coef : {0:.2f}'.format(blood_tissue_h2o_part_coef), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + rf_block_duration : {0:.4f}'.format(rf_block_duration), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + tag_delay_ms : {0:d}'.format(tag_delay_ms), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + min_tr_ms : {0:d}'.format(min_tr_ms), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + num_rf_blocks : {0:d}'.format(num_rf_blocks), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + num_slices : {0:d}'.format(num_slices), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + tagging_time : {0:d}'.format(tagging_time), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + slice_time : {0:.8f}'.format(slice_time), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + r1a : {0:.8f}'.format(r1a), 'red', 'on_grey', attrs=['bold'] )
        cprint('  + tagging_efficiency : {0:.2f}'.format(tagging_efficiency), 'red', 'on_grey', attrs=['bold'] )
        if i_file_m0 is not None:
            cprint('  + using PD image for m0 normalization', 'red', 'on_grey', attrs=['bold'] )
        else:
            cprint('  + using control image normalization', 'red', 'on_grey', attrs=['bold'] )
        cprint("=======================================================", 'red', 'on_grey', attrs=['bold'] )




        #    we are creating a slice timing adjustment for the total delay
        #    which includes the post-tagging delay plus the timing of each slice acquistion from bottom to top of the brain
        #    don't use with slice timing corrected input volumes [EITHER/OR]. This method used for ASLToolbox comparison
        #    we might consider using AFNI's slice timing correction instead and removing this slice_time_array
        #
        #    create a slicetimearray of ones (64*64,24)
        #    also add delay time aka post_label delay
        slice_time_array = np.ones(asl_xdim*asl_ydim*asl_zdim,dtype=float).reshape(asl_xdim,asl_ydim,asl_zdim)
        slice_weight = (np.arange(0,asl_zdim,1)*slice_time + post_tag_delay_time)
        for slice_indx in range(asl_zdim):
            slice_time_array[:,:,slice_indx] = slice_weight[slice_indx]



        # need mean images to store in below
        mean_cbf_vol = np.zeros(asl_xdim*asl_ydim*asl_zdim,dtype=float).reshape(asl_xdim,asl_ydim,asl_zdim)
        # array for all cbf volumes
        #all_cbf_vol = np.zeros(asl_xdim*asl_ydim*asl_zdim*img_pairs,dtype=float).reshape(asl_xdim,asl_ydim,asl_zdim,img_pairs)



        ##################################################
        # conduct CBF estimates...
        ##################################################
        print("Estimating CBF...")
        # foreach label-control pair:
        # for now we'll assume that input is label first, control second
        # this per pair method is for comparison to ASLToolbox which implements it this way.
        # Consider averaging control and labelled images separately and calculating cbf once on those
        for pair_indx in np.arange(0,asl_tdim,2):

            #	subtract control - label and call perfusion image
            pair_perf = asl_AfniImg.data[:,:,:,pair_indx+1] - asl_AfniImg.data[:,:,:,pair_indx]

            #	Divide perfusion image by control image(M0)
            # since i loaded asl data as int16 to save memory i use true divide here so that result is float
            # m0_vol is masked array so only relevant data is used in division
            if i_file_m0 is None:
                pair_perf_norm = np.true_divide(pair_perf,asl_AfniImg.data[:,:,:,pair_indx+1])
            else:
                pair_perf_norm = np.true_divide(pair_perf,m0_vol)

            #	calculate efftcbf with formula for that pair "efftcbf=6000*1000*lambda*efftcbf*r1a./(2*labeff* (exp(-omega*r1a)-exp( -1*(Labeltime+omega)*r1a) ) );"
            # this formula format is to match ASLToolbox
            # [ml/100g/min]
            cbf = 6000*1000*blood_tissue_h2o_part_coef*pair_perf_norm*r1a/(2*tagging_efficiency*(np.exp(-slice_time_array*r1a)-np.exp(-1*(tagging_time+slice_time_array)*r1a)))

            #	add this to mean cbf image
            mean_cbf_vol = mean_cbf_vol + cbf

            # all cbf volumes
            #all_cbf_vol[:,:,:,pair_indx/2] = cbf



        # finally, write out mean cbf volume
        # TODO: if i_file_m0 not used the cbf data doesn't get masked.
        print("Writing Mean CBF Volume...")

        if o_path is None:
            o_path = i_path

        # remove nan's
        mean_cbf_vol = np.nan_to_num(mean_cbf_vol)
        # get mean cbf and reshape
        mean_cbf_vol /= img_pairs
        # create output image from input image
        if i_file_m0 is not None:
            cbf_img = utils.AfniImage(np.array(mean_cbf_vol,dtype=np.float32), m0_AfniImg.header)
        else:
            cbf_img = utils.AfniImage(np.array(mean_cbf_vol,dtype=np.float32), asl_AfniImg.header)
        # save the image to file
        cbf_img.save(os.path.join(o_path,o_prfx + ".mean.cbf"), gz=True)

        #        nibabel NIFTI
        #        o_img = nb.Nifti1Image(np.array(mean_cbf_vol,dtype=np.float32), asl_niftiImg.affine, asl_niftiImg.header)
        #        # clear the AFNI extension from this output image
        #        o_img.header.extensions = o_img.header.exts_klass()
        #        # save the image to file
        #        nb.save(o_img, os.path.join(work_path,o_prfx + ".mean.cbf.nii.gz"))



    except:

        raise





def qa_spatial_normalization(subj, study):

    """
    QA display of pcasl CBF spatial normalization for inspection.

    This display CBF over subject's non-linear structural warp.

    :param subj: utils.Subj
    :param study: utils.Study
    """

    import shutil, glob, os
    from mayerlab import utils

    workspace_path = "/tmp"

    #    try doing this in machine tmp since this is all temporary output
    os.chdir(workspace_path)

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : pcasl : spatial normalization check" )

    # check
    if not os.path.isfile(subj.get_path(study,"anat") + study.label + "." + subj.ursi.short+".anat.sksp.qwarp+tlrc.HEAD"):
        raise utils.CmdProcError("Missing T1 qwarp")
    
    if not os.path.isfile(os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.mean.cbf.fwhm6+tlrc.HEAD") ):
        raise utils.CmdProcError("Missing CBF")


    # get either TT_N27 or qwarp
    try:
        for filename in glob.glob(os.path.join(subj.get_path(study,"anat"), study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc*")):
            shutil.copy(filename, workspace_path)
    except:
        raise    
        
    # get already normalized asl
    try:
        for filename in glob.glob( os.path.join( subj.get_path(study), 'pcasl', study.label + "." + subj.ursi.short+".pCASL.mean.cbf.fwhm6+tlrc.*") ):
            shutil.copy(filename, workspace_path)
    except:
        raise


    cmd = utils.Cmd("3dcalc -overwrite -a " + subj.get_path(study) + "pcasl/" + study.label + "." + subj.ursi.short + ".pCASL.3dregToPD+orig'[0]' -expr 'a' -prefix /tmp/VERIFY.pcasl" )
    cmd.run()
    
    cmd = utils.Cmd("3dAllineate " +
            " -master " + subj.get_path(study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.sksp+tlrc " +
            " -1Dmatrix_apply " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + ".pCASL_epi2tlrc.aff12.1D" +
            " -input /tmp/VERIFY.pcasl+orig " +
            " -floatize" +
            " -interp trilinear" +
            " -final wsinc5"+
            " -prefix /tmp/VERIFY.pcasl.tmp" )

    cmd.run()
    
    cmd = utils.Cmd("3dedge3 -overwrite -prefix /tmp/VERIFY.pcasl.tmp.edge -input /tmp/VERIFY.pcasl.tmp+tlrc.HEAD" )
    cmd.run()


    # scale the volume for better contrast
    cmd = utils.Cmd("@ScaleVolume -input /tmp/VERIFY.pcasl.tmp.edge+tlrc -prefix /tmp/VERIFY.pcasl.tmp.edge.sc -perc_clip 2 98 -val_clip 6 249")
    cmd.run()

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : pcasl : spatial normalization check : Base" )
    cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subj.ursi.short + '.anat.sksp.qwarp+tlrc" -com "SWITCH_OVERLAY VERIFY.pcasl.tmp.edge.sc" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ')
    cmd.run()
    
    
    
    
    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : pcasl : spatial normalization check : CBF" )
    # mask out using underlay
    utils.Cmd("3dcalc -overwrite -a "+ study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc -b " + study.label + "." + subj.ursi.short+".pCASL.mean.cbf.fwhm6+tlrc -expr 'step(a)*b' -prefix VERIFY.pcasl.cbf.tmp").run()
    # open afni
    utils.Cmd('afni -no_detach -q -noplugins -posfunc -com "SWITCH_UNDERLAY ' + study.label + '.' + subj.ursi.short + '.anat.sksp.qwarp+tlrc" -com "SWITCH_OVERLAY VERIFY.pcasl.cbf.tmp+tlrc" -com "SEE_OVERLAY +" -com "SET_THRESHNEW A 9.999" -com "SET_FUNC_RANGE A 80" -com "SET_SPM_XYZ 6 12 14" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_PBAR_NUMBER A.14" ').run()
    
    
    #    clean up
    utils.purge_paths("/tmp/" + study.label + "." +  subj.ursi.short + ".*", "/tmp/*pcasl*")




def qa_mask(subj, study):

    from mayerlab import utils
    import os

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : pcasl : mask check" )

    os.chdir(subj.get_path(study) + "pcasl")
    
    utils.Cmd('3dAutomask -overwrite -q -clfrac 0.70 -peels 2 -nbhrs 16 -apply_prefix ' + study.label + '.' + subj.ursi.short + '.PD.newmask+orig.HEAD ' + study.label + '.' + subj.ursi.short + '.PD+orig.HEAD').run()
    
    cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subj.ursi.short + '.PD+orig.HEAD" -com "SWITCH_OVERLAY ' + study.label + '.' + subj.ursi.short + '.PD.newmask+orig" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ')
    cmd.run()

    utils.purge_paths("*newmask*")
    
    

    

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

