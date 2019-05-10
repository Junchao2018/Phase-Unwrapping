# Phase-Unwrapping
Phase unwrapping using deep learning

Please cite the following paper and patent when you use it:
Paper:
Junchao Zhang, Xiaobo Tian, Jianbo Shao, Haibo Luo, and Rongguang Liang.
"Phase unwrapping in optical metrology via denoised and convolutional
segmentation networks,"  Optics Express 27(10), 14903-14912 (2019).

Patent:
Rongguang Liang, Junchao Zhang, Xiaobo Tian, and Jianbo Shao. Phase
unwrapping using segmentation. (U.S. Provisional Patent Application, No.
62/768,624)

-----------------------------------------------------------------------------------------------
The segmentation network can be used for phase discontinuity extraction as well. For more detail,
please read our paper.
-----------------------------------------------------------------------------------------------
Denoised Network: 
it is based on U-Net with removing pooling and upsampling layers. Besides, two observations
are used to train the network.
-----------------------------------------------------------------------------------------------
