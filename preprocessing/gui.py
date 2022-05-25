#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:21:40 2021

@author: josef
"""

import sys

from PyQt5.QtWidgets import (QCheckBox, QLabel, QDialog, QDialogButtonBox,  
                             QHBoxLayout, QVBoxLayout, QApplication)

        
    


class UserConfirm(QDialog):     

    def __init__(self, title, msg):
        super(UserConfirm, self).__init__()

        self.title = title
        self.msg = msg
        
        self.initUI()
        
        
    def initUI(self):

        vbox_selections = QVBoxLayout()      
        
        msg = QLabel(self.msg)
        vbox_selections.addWidget(msg)
        vbox_selections.addSpacing(10)
        
        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(vbox_selections)
        
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        buttonBox.centerButtons()
        hbox_buttons = QHBoxLayout()
        hbox_buttons.addStretch(1)
        hbox_buttons.addWidget(buttonBox)

        vbox.addSpacing(20)
        vbox.addLayout(hbox_buttons)
        
        
        self.setLayout(vbox)

        self.setGeometry(300, 300, 300, 150)
        self.setWindowTitle(self.title)
        self.show()
        
        
    def accept(self):            
        super(UserConfirm, self).accept()

    def reject(self):        
        super(UserConfirm, self).reject()




def get_user_confirm(title, msg):
    
    try:
        _ = QApplication(sys.argv)
        gui = UserConfirm(title, msg)
        gui.show()
        
    except Exception as e:
        print(e)

    finally:
        result = gui.exec_()
        if result == QDialog.Accepted:
            return True
        else:
            return False
        
        
    





class UserSelectMultiple(QDialog):     

    def __init__(self, title, msg, inputList, preselectedList):
        #super().__init__()
        super(UserSelectMultiple, self).__init__()

        self.title = title
        self.msg = msg
        self.inputList = inputList
        self.preselectedList = preselectedList
        
        # create checkboxes in loop. Other forms of dynamic initialization did not allow retrieval of proper state
        self.cbs = []
        for i, input in enumerate(self.inputList):
            self.cbs.append(QCheckBox(str(i)))
        
        
        self.initUI()
        
        
    def initUI(self):

        vbox_selections = QVBoxLayout()      
        
        msg = QLabel(self.msg)
        vbox_selections.addWidget(msg)
        vbox_selections.addSpacing(10)
        
        
        for i, input in enumerate(self.inputList):
            self.cbs[i].setText(input)
            if i in self.preselectedList:
                self.cbs[i].toggle()
            vbox_selections.addWidget(self.cbs[i])


        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(vbox_selections)
        
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        buttonBox.centerButtons()
        hbox_buttons = QHBoxLayout()
        hbox_buttons.addStretch(1)
        hbox_buttons.addWidget(buttonBox)

        vbox.addSpacing(20)
        vbox.addLayout(hbox_buttons)
        
        
        self.setLayout(vbox)

        self.setGeometry(300, 300, 300, 150)
        self.setWindowTitle(self.title)
        self.show()
        
        
    def accept(self):            
        super(UserSelectMultiple, self).accept()

    def reject(self):        
        super(UserSelectMultiple, self).reject()
        
    def get_output(self):
        return_list = []
        for i, cb in enumerate(self.cbs):
            if cb.checkState():
                #return_list.append(cb.text())
                return_list.append(i)
                
        return return_list
    



def get_user_selections(title, msg, inputList, preselectedList=[]):
    
    try:
        _ = QApplication(sys.argv)
        gui = UserSelectMultiple(title, msg, inputList, preselectedList)
        gui.show()
        
    except Exception as e:
        print(e)

    finally:
        result = gui.exec_()
        if result == QDialog.Accepted:
            return gui.get_output()
        else:
            return []
        
        
    