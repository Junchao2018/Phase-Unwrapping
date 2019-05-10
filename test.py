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
import tensorflow as tf
import model as model
import numpy as np
import h5py
import math
import os
import scipy.io

MODEL_SAVE_PATH = '/xdisk/junchaozhang/model/'
CH = 1
IMG_SIZE = (275,275)



def test(test_data,ground_truth):
    with tf.Graph().as_default():
        x = tf.placeholder(tf.float32,[None,
                                       IMG_SIZE[0],
                                       IMG_SIZE[1],
                                       CH])

        logits = model.forward(x, False,with_dropout=False,keep_prob=0.5,batch_size=1)
        logits_reshape = tf.reshape(logits, [-1, model.Num_Classes])
        prob = tf.nn.softmax(logits_reshape, -1)
        prediction = tf.argmax(prob, -1)
        y = tf.reshape(prediction, IMG_SIZE)


        saver = tf.train.Saver()
        
        with tf.Session() as sess:
            ckpt = tf.train.get_checkpoint_state(MODEL_SAVE_PATH)
            if ckpt and ckpt.model_checkpoint_path:
                    saver.restore(sess,ckpt.model_checkpoint_path)
                    global_step = ckpt.model_checkpoint_path.split('/')[-1].split('-')[-1]
                    print(global_step)


                    Num = test_data.shape[0]
                    ImgOut = np.zeros([Num,IMG_SIZE[0],IMG_SIZE[1]],dtype = np.int64)
                    for i in range(Num):
                        img = sess.run(y, feed_dict={x: test_data[i:i+1,:,:,:]})
                        ImgOut[i,:,:] = np.array(img)

                    scipy.io.savemat('results.mat', {'mydata': ImgOut})
            else:
                print("No checkpoint is found.")
                return

if __name__=='__main__':
    data = h5py.File('/xdisk/junchaozhang/Data/imtest.mat')
    input_data = data["test_inputs"]
    output_data = data["test_outputs"]
    input_npy = np.transpose(input_data)
    output_npy = np.transpose(output_data)

    #######################
    input_npy = np.reshape(input_npy, (input_npy.shape[0], input_npy.shape[1], input_npy.shape[2], 1))
    output_npy = np.reshape(output_npy,(output_npy.shape[0], output_npy.shape[1], output_npy.shape[2],1))
    print(np.shape(input_npy))
    print(np.shape(output_npy))

    test(input_npy,output_npy)
        
        

