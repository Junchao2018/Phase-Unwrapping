# -*- coding: utf-8 -*-
"""
========================================================================
This network is  for Phase Unwrapping, Version 1.0
Copyright(c) 2019  Junchao Zhang, Xiaobo Tian, Jianbo Shao, Haibo Luo,
and Rongguang Liang
All Rights Reserved.
----------------------------------------------------------------------
Permission to use, copy, or modify this software and its documentation
for educational and research purposes only and without fee is here
granted, provided that this copyright notice and the original authors'
names appear on all copies and supporting documentation. This program
shall not be used, rewritten, or adapted as the basis of a commercial
software or hardware product without first obtaining permission of the
authors. The authors make no representations about the suitability of
this software for any purpose. It is provided "as is" without express
or implied warranty.
----------------------------------------------------------------------
Please cite the following paper and patent when you use it:
Paper:
Junchao Zhang, Xiaobo Tian, Jianbo Shao, Haibo Luo, and Rongguang Liang.
"Phase unwrapping in optical metrology via denoised and convolutional
segmentation networks,"  Optics Express 27(10), 14903-14912 (2019).

Patent:
Rongguang Liang, Junchao Zhang, Xiaobo Tian, and Jianbo Shao. Phase
unwrapping using segmentation. (U.S. Provisional Patent Application, No.
62/768,624)
========================================================================
"""
import numpy as np

def data_augmentation(imagein, mode):
    image = np.transpose(imagein,(1,2,3,0))
    if mode == 0:
        # original
        return imagein
    elif mode == 1:
        # flip up and down
        image = np.flipud(image)
    elif mode == 2:
        # rotate counterwise 90 degree
        image = np.rot90(image)
    elif mode == 3:
        # rotate 90 degree and flip up and down
        image = np.rot90(image)
        image = np.flipud(image)
    elif mode == 4:
        # rotate 180 degree
        image = np.rot90(image, k=2)
    elif mode == 5:
        # rotate 180 degree and flip
        image = np.rot90(image, k=2)
        image = np.flipud(image)
    elif mode == 6:
        # rotate 270 degree
        image = np.rot90(image, k=3)
    elif mode == 7:
        # rotate 270 degree and flip
        image = np.rot90(image, k=3)
        image = np.flipud(image)
    imageout = np.transpose(image,(3,0,1,2))
    return imageout