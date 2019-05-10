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
import math

Num_Classes = 15


def max_pool(inputs, name):
    with tf.variable_scope(name) as scope:
        value, index = tf.nn.max_pool_with_argmax(tf.to_double(inputs), ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1],
                                                  padding='SAME', name=scope.name)
    return tf.to_float(value), index, inputs.get_shape().as_list()


def initialization(k, c):
    std = math.sqrt(2. / (k ** 2 * c))
    return tf.truncated_normal_initializer(stddev=std)


def conv_layer(bottom, name, shape, is_training):
    with tf.variable_scope(name) as scope:
        w = tf.get_variable('weights', shape, initializer=initialization(shape[0], shape[2]))
        conv = tf.nn.conv2d(bottom, w, [1,1,1,1], padding='SAME')
        b = tf.get_variable('bias', shape[3], initializer=tf.constant_initializer(0.0))
        bias = tf.nn.bias_add(conv,b)
        bias = tf.layers.batch_normalization(bias, training=is_training)
        conv_out = tf.nn.relu(bias)
    return conv_out



def up_sampling(pool, ind, output_shape, batch_size, name=None):
    with tf.variable_scope(name):
        pool_ = tf.reshape(pool, [-1])
        batch_range = tf.reshape(tf.range(batch_size, dtype=ind.dtype), [tf.shape(pool)[0], 1, 1, 1])
        b = tf.ones_like(ind) * batch_range
        b = tf.reshape(b, [-1, 1])
        ind_ = tf.reshape(ind, [-1, 1])
        ind_ = tf.concat([b, ind_], 1)
        ret = tf.scatter_nd(ind_, pool_, shape=[batch_size, output_shape[1] * output_shape[2] * output_shape[3]])
        ret = tf.reshape(ret, [tf.shape(pool)[0], output_shape[1], output_shape[2], output_shape[3]])
        return ret


def forward(inputs,is_training,with_dropout,keep_prob, batch_size, bayes = True):
    with tf.variable_scope('forward'):
        conv1_1 = conv_layer(inputs,'conv1_1',[3, 3, 1, 64], is_training)
        conv1_2 = conv_layer(conv1_1,'conv1_2',[3, 3, 64, 64], is_training)
        pool1, pool1_index, shape_1 = max_pool(conv1_2, 'pool1')

        conv2_1 = conv_layer(pool1,'conv2_1',[3, 3, 64, 128], is_training)
        conv2_2 = conv_layer(conv2_1,'conv2_2',[3, 3, 128, 128], is_training)
        pool2, pool2_index, shape_2 = max_pool(conv2_2, 'pool2')

        conv3_1 = conv_layer(pool2,'conv3_1',[3, 3, 128, 256], is_training)
        conv3_2 = conv_layer(conv3_1,'conv3_2',[3, 3, 256, 256], is_training)
        conv3_3 = conv_layer(conv3_2,'conv3_3',[3, 3, 256, 256], is_training)
        pool3, pool3_index, shape_3 = max_pool(conv3_3, 'pool3')

        if bayes:
            dropout1 = tf.layers.dropout(pool3, rate = (1-keep_prob), training = with_dropout, name='dropout1')
            conv4_1 = conv_layer(dropout1, 'conv4_1', [3, 3, 256, 512], is_training)
        else:
            conv4_1 = conv_layer(pool3, 'conv4_1', [3, 3, 256, 512], is_training)
        conv4_2 = conv_layer(conv4_1, 'conv4_2', [3, 3, 512, 512], is_training)
        conv4_3 = conv_layer(conv4_2, 'conv4_3', [3, 3, 512, 512], is_training)
        pool4, pool4_index, shape_4 = max_pool(conv4_3, 'pool4')

        if bayes:
            dropout2 = tf.layers.dropout(pool4, rate = (1-keep_prob), training = with_dropout, name='dropout2')
            conv5_1 = conv_layer(dropout2, 'conv5_1', [3, 3, 512, 512], is_training)
        else:
            conv5_1 = conv_layer(pool4, 'conv5_1', [3, 3, 512, 512], is_training)
        conv5_2 = conv_layer(conv5_1, 'conv5_2', [3, 3, 512, 512], is_training)
        conv5_3 = conv_layer(conv5_2, 'conv5_3', [3, 3, 512, 512], is_training)
        pool5, pool5_index, shape_5 = max_pool(conv5_3, 'pool5')
        
        #Decoder Process
        if bayes:
            dropout3 = tf.layers.dropout(pool5, rate = (1-keep_prob), training = with_dropout, name='dropout3')
            deconv5_1 = up_sampling(dropout3, pool5_index, shape_5, batch_size, name = 'unpool_5')
        else:
            deconv5_1 = up_sampling(pool5, pool5_index, shape_5, batch_size, name = 'unpool_5')
        deconv5_2 = conv_layer(deconv5_1, 'deconv5_2', [3, 3, 512, 512], is_training)
        deconv5_3 = conv_layer(deconv5_2, 'deconv5_3', [3, 3, 512, 512], is_training)
        deconv5_4 = conv_layer(deconv5_3, 'deconv5_4', [3, 3, 512, 512], is_training)

        if bayes:
            dropout4 = tf.layers.dropout(deconv5_4, rate = (1-keep_prob), training = with_dropout, name='dropout4')
            deconv4_1 = up_sampling(dropout4, pool4_index, shape_4, batch_size, name = 'unpool_4')
        else:
            deconv4_1 = up_sampling(deconv5_4, pool4_index, shape_4, batch_size, name = 'unpool_4')
        deconv4_2 = conv_layer(deconv4_1, 'deconv4_2', [3, 3, 512, 512], is_training)
        deconv4_3 = conv_layer(deconv4_2, 'deconv4_3', [3, 3, 512, 512], is_training)
        deconv4_4 = conv_layer(deconv4_3, 'deconv4_4', [3, 3, 512, 256], is_training)


        if bayes:
            dropout5 = tf.layers.dropout(deconv4_4, rate = (1-keep_prob), training = with_dropout, name='dropout5')
            deconv3_1 = up_sampling(dropout5, pool3_index, shape_3, batch_size, name = 'unpool_3')
        else:
            deconv3_1 = up_sampling(deconv4_4, pool3_index, shape_3, batch_size, name = 'unpool_3')
        deconv3_2 = conv_layer(deconv3_1, 'deconv3_2', [3, 3, 256, 256], is_training)
        deconv3_3 = conv_layer(deconv3_2, 'deconv3_3', [3, 3, 256, 256], is_training)
        deconv3_4 = conv_layer(deconv3_3, 'deconv3_4', [3, 3, 256, 128], is_training)

        if bayes:
            dropout6 = tf.layers.dropout(deconv3_4, rate = (1-keep_prob), training = with_dropout, name='dropout6')
            deconv2_1 = up_sampling(dropout6, pool2_index, shape_2, batch_size, name = 'unpool_2')
        else:
            deconv2_1 = up_sampling(deconv3_4, pool2_index, shape_2, batch_size, name = 'unpool_2')
        deconv2_2 = conv_layer(deconv2_1, 'deconv2_2', [3, 3, 128, 128], is_training)
        deconv2_3 = conv_layer(deconv2_2, 'deconv2_3', [3, 3, 128, 64], is_training)

        deconv1_1 = up_sampling(deconv2_3, pool1_index, shape_1, batch_size, name = 'unpool_1')
        deconv1_2 = conv_layer(deconv1_1, 'deconv1_2', [3, 3, 64, 64], is_training)
        deconv1_3 = conv_layer(deconv1_2, 'deconv1_3', [3, 3, 64, 64], is_training)


        with tf.variable_scope('conv_classifier') as scope:
            kernel = tf.get_variable('weights', [1,1,64,Num_Classes], initializer=initialization(1, 64))
            conv = tf.nn.conv2d(deconv1_3, kernel, [1,1,1,1], padding='VALID')
            bias = tf.get_variable('bias', Num_Classes, initializer=tf.constant_initializer(0.0))
            logits = tf.nn.bias_add(conv,bias)
    return logits