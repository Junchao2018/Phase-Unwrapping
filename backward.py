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
import model
import os
import numpy as np
import data_augmentation as DA
import h5py



BATCH_SIZE = 3
IMG_SIZE = (275,275)
IMG_CHANNEL = 1
MAX_EPOCH = 1000


LEARNING_RATE_BASE = 0.0001
LEARNING_RATE_DECAY = 0.99
MODEL_SAVE_PATH = './model/'
MODEL_NAME = 'Phase_Model'




def normal_loss(logits, labels, number_class):
    label_flatten = tf.to_int64(tf.reshape(labels, [-1]))
    logits_reshape = tf.reshape(logits, [-1, number_class])

    Num_Label = tf.bincount(tf.to_int32(label_flatten),minlength = number_class)
    #
    prob = tf.nn.softmax(logits_reshape, -1)
    prediction = tf.argmax(prob, -1)
    Num_prediction = tf.bincount(tf.to_int32(prediction),minlength = number_class)

    Num_diff = tf.abs(Num_Label - Num_prediction)
    frequency = (tf.add(tf.to_float(Num_diff),0.0001))/tf.reduce_sum(tf.add(tf.to_float(Num_diff),0.0001))

    label_onehot = tf.one_hot(label_flatten, depth=number_class)
    cross_entropy = tf.nn.weighted_cross_entropy_with_logits(targets=label_onehot, logits=logits_reshape,
                                                             pos_weight=frequency)

    #cross_entropy = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=label_flatten, logits=logits_reshape,
    #                                                               name='normal_cross_entropy')
    cross_entropy_mean = tf.reduce_mean(cross_entropy, name='cross_entropy')
    correct_prediction = tf.equal(tf.argmax(logits_reshape, -1), label_flatten)
    accuracy = tf.reduce_mean(tf.to_float(correct_prediction))

    return cross_entropy_mean, accuracy



def backward(inputs,labels,train_num):
    with tf.Graph().as_default() as g:
        with tf.name_scope('input'):
            x = tf.placeholder(tf.float32, [None, IMG_SIZE[0],IMG_SIZE[1],IMG_CHANNEL])
            y_ = tf.placeholder(tf.int64, [None, IMG_SIZE[0],IMG_SIZE[1], 1])
            with_dropout = tf.placeholder(tf.bool, name="with_dropout")
            keep_prob = tf.placeholder(tf.float32, shape=None, name="keep_rate")
            batch_size = tf.placeholder(tf.int64, shape=[], name="batch_size")

        global_step = tf.Variable(0,trainable=False)
        learning_rate = tf.train.exponential_decay(LEARNING_RATE_BASE,global_step,
                                                   train_num//BATCH_SIZE,
                                                   LEARNING_RATE_DECAY, staircase=True)


        logits = model.forward(x,True,with_dropout,keep_prob,batch_size)
        with tf.name_scope('loss'):
            loss, accuracy = normal_loss(logits,y_,model.Num_Classes)

        
        optimizer = tf.train.AdamOptimizer(learning_rate)
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(update_ops):
            train_op = optimizer.minimize(loss,global_step=global_step)


        # Save model
        saver = tf.train.Saver(max_to_keep = 50)
        epoch = 0
        with tf.Session() as sess:
            tf.global_variables_initializer().run()
            ckpt = tf.train.get_checkpoint_state(MODEL_SAVE_PATH)
            if ckpt and ckpt.model_checkpoint_path:
                saver.restore(sess,ckpt.model_checkpoint_path)
                epoch = int(ckpt.model_checkpoint_path.split('/')[-1].split('_')[-1].split('-')[-2])
            
            while epoch<MAX_EPOCH:
                max_step = train_num//BATCH_SIZE
                listtmp = np.random.permutation(train_num)
                j = 0
                for i in range(max_step):
                    file =open("loss.txt",'a')
                    ind = listtmp[j:j+BATCH_SIZE]
                    j = j + BATCH_SIZE
                    xs = inputs[ind,:,:,:]
                    ys = labels[ind,:,:,:]
                    mode = np.random.permutation(8)
                    xs = DA.data_augmentation(xs,mode[0])
                    ys = DA.data_augmentation(ys,mode[0])
                    _,loss_v,accu, step = sess.run([train_op,loss,accuracy,global_step],feed_dict={x: xs, y_: ys,
                                                                                                   with_dropout: True,
                                                                                                   keep_prob: 0.5,
                                                                                                   batch_size: BATCH_SIZE})
                    file.write("Epoch: %d, After [%d / %d] training, the batch loss is %g, the accuracy is %g.\n" %(epoch,i+1,max_step,loss_v,accu))
                    file.close()
                saver.save(sess,os.path.join(MODEL_SAVE_PATH,MODEL_NAME+'_epoch_'+str(epoch+1)),global_step = global_step)
                epoch +=1

if __name__=='__main__':
    data = h5py.File('/xdisk/junchaozhang/Data/imtrain.mat')
    input_data = data["training_inputs"]
    output_data = data["training_outputs"]
    input_npy = np.transpose(input_data)
    output_npy = np.transpose(output_data)
    input_npy = np.reshape(input_npy, (input_npy.shape[0], input_npy.shape[1], input_npy.shape[2], 1))
    output_npy = np.reshape(output_npy,(output_npy.shape[0], output_npy.shape[1], output_npy.shape[2],1))

    print(np.shape(input_npy))
    print(np.shape(output_npy))

    train_num = input_npy.shape[0]
    backward(input_npy,output_npy,train_num)
