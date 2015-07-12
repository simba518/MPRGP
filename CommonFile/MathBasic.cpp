#include "MathBasic.h"

USE_PRJ_NAMESPACE

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

sizeType ScalarUtil<sizeType>::scalar_max=numeric_limits<sizeType>::max();
sizeType ScalarUtil<sizeType>::scalar_eps=1;

OmpSettings OmpSettings::_ompSettings;

const scalarD ErfCalc::cof[28] = 
{
    -1.3026537197817094, 6.4196979235649026e-1,
    1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
    3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
    -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
    6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
    9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
    -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17
};
