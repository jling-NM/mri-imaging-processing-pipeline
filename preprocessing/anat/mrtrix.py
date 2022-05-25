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



def degibbs(subj):
    """
    Remove gibbs ringing from T1w
    """
    anat_path = subj.derivatives_path(bidsType='anat')
    
    i_file = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    o_file = os.path.join(anat_path, subj.bids_file_join('T1w.dgibbs.nii.gz'))
    
    if not os.path.isfile(o_file):
        utils.Cmd("mrdegibbs " + i_file + " " + o_file).echo()
