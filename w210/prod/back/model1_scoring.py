#!/usr/bin/env python
import calendar
from datetime import datetime
from collections import namedtuple
import re
import sys
import time
import os

import numpy as np
import pandas as pd

PY3 = sys.version_info[0] == 3
if PY3:
    string_types = str,
    text_type = str
    long_type = int
else:
    string_types = basestring,
    text_type = unicode
    long_type = long

def predict(row):
    round_VEGFR = np.float32(row[u'VEGFR'])
    round_drug_id__308 = np.float32(row[u'drug_id__308'])
    round_387266 = np.float32(row[u'387266'])
    round_cna_princ10 = np.float32(row[u'cna_princ10'])
    round_Num_Targets = np.float32(row[u'Num_Targets'])
    round_CSF1R = np.float32(row[u'CSF1R'])
    round_cna_princ32 = np.float32(row[u'cna_princ32'])
    round_SRC = np.float32(row[u'SRC'])
    round_cna_princ40 = np.float32(row[u'cna_princ40'])
    round_cna_princ34 = np.float32(row[u'cna_princ34'])
    round_ALK = np.float32(row[u'ALK'])
    round_cna_princ18 = np.float32(row[u'cna_princ18'])
    round_cna_princ19 = np.float32(row[u'cna_princ19'])
    round_cna_princ31 = np.float32(row[u'cna_princ31'])
    round_cna_princ14 = np.float32(row[u'cna_princ14'])
    round_cna_princ15 = np.float32(row[u'cna_princ15'])
    round_cna_princ17 = np.float32(row[u'cna_princ17'])
    round_drug_id__35 = np.float32(row[u'drug_id__35'])
    round_cna_princ38 = np.float32(row[u'cna_princ38'])
    round_cna_princ13 = np.float32(row[u'cna_princ13'])
    round_drug_id__1029 = np.float32(row[u'drug_id__1029'])
    round_drug_id__1023 = np.float32(row[u'drug_id__1023'])
    round_2019 = np.float32(row[u'2019'])
    round_drug_id__1021 = np.float32(row[u'drug_id__1021'])
    round_RET = np.float32(row[u'RET'])
    round_347733 = np.float32(row[u'347733'])
    round_None = np.float32(row[u'None'])
    round_drug_id__193 = np.float32(row[u'drug_id__193'])
    round_ARID1A_mut = np.float32(row[u'ARID1A_mut'])
    round_FLT3 = np.float32(row[u'FLT3'])
    round_AFF4_mut = np.float32(row[u'AFF4_mut'])
    round_drug_id__312 = np.float32(row[u'drug_id__312'])
    round_4633 = np.float32(row[u'4633'])
    round_drug_id__199 = np.float32(row[u'drug_id__199'])
    round_drug_id__6 = np.float32(row[u'drug_id__6'])
    round_MSI_FACTOR = np.float32(row[u'MSI_FACTOR'])
    round_KIT = np.float32(row[u'KIT'])
    round_FGFR1 = np.float32(row[u'FGFR1'])
    round_TIE2 = np.float32(row[u'TIE2'])
    round_1063 = np.float32(row[u'1063'])
    round_drug_id__254 = np.float32(row[u'drug_id__254'])
    round_drug_id__281 = np.float32(row[u'drug_id__281'])
    round_cna_princ33 = np.float32(row[u'cna_princ33'])
    round_cna_princ25 = np.float32(row[u'cna_princ25'])
    round_cna_princ27 = np.float32(row[u'cna_princ27'])
    round_cna_princ20 = np.float32(row[u'cna_princ20'])
    round_cna_princ23 = np.float32(row[u'cna_princ23'])
    round_cna_princ22 = np.float32(row[u'cna_princ22'])
    round_cna_princ2 = np.float32(row[u'cna_princ2'])
    round_drug_id__1194 = np.float32(row[u'drug_id__1194'])
    round_cna_princ29 = np.float32(row[u'cna_princ29'])
    round_cna_princ28 = np.float32(row[u'cna_princ28'])
    round_cna_princ1 = np.float32(row[u'cna_princ1'])
    round_VEGFR1 = np.float32(row[u'VEGFR1'])
    round_VEGFR2 = np.float32(row[u'VEGFR2'])
    round_ALK5 = np.float32(row[u'ALK5'])
    round_cna_princ8 = np.float32(row[u'cna_princ8'])
    round_drug_id__30 = np.float32(row[u'drug_id__30'])
    round_cna_princ6 = np.float32(row[u'cna_princ6'])
    round_cna_princ7 = np.float32(row[u'cna_princ7'])
    round_cna_princ4 = np.float32(row[u'cna_princ4'])
    round_cna_princ5 = np.float32(row[u'cna_princ5'])
    round_drug_id__298 = np.float32(row[u'drug_id__298'])
    round_cna_princ3 = np.float32(row[u'cna_princ3'])
    round_TP53_mut = np.float32(row[u'TP53_mut'])
    round_cna_princ9 = np.float32(row[u'cna_princ9'])
    return sum([
         2.6200358,
         0.00023022528496268156997 * (round_cna_princ15),
         0.00053821001342795472388 * (round_cna_princ19),
        -0.00041471436106392029928 * (round_cna_princ31),
         -0.0015668989121138876194 * (round_cna_princ38),
           -0.12366105299613171231 * (round_347733),
           -0.31257858115864628479 * (round_4633),
           -0.11944903885526887699 * (round_2019),
           0.091864431015139622061 * (round_drug_id__6),
          -0.012625470495153976597 * (round_drug_id__30),
           -0.37903175747727058109 * (round_drug_id__35),
            0.61123394044325141117 * (round_drug_id__193),
           -0.40532689945545247179 * (round_drug_id__308),
           -0.97803994995082477182 * (round_drug_id__312),
            -0.1242826606927681099 * (round_MSI_FACTOR),
           -0.20854423789886750162 * (round_AFF4_mut),
           -0.03143723810763280474 * (round_TP53_mut),
          -0.037395196938429596523 * (round_SRC),
           -0.35110369480625480687 * (round_TIE2),
          0.0091575937723347693858 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5 and 
                                     round_TIE2 <= 0.5),
           0.096733252853388768111 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_None > 2.0 and 
                                     round_VEGFR2 <= 0.5),
          -0.014397353659776528437 * (round_drug_id__193 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5 and 
                                     round_None > 2.5),
           -0.44124464731168078213 * (round_None <= 1.5),
           0.045932022021736268469 * (round_cna_princ6 > 0.28661248087882996 and 
                                     round_cna_princ7 > -23.252197265625 and 
                                     round_1063 <= 0.5 and 
                                     round_FGFR1 <= 0.5),
            0.38289567744623437262 * (round_cna_princ34 <= -0.14065495133399963 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_Num_Targets <= 7.0 and 
                                     round_VEGFR2 <= 0.5),
            0.40577335736104114261 * (round_cna_princ34 > -0.14065495133399963 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_Num_Targets <= 7.0 and 
                                     round_VEGFR2 <= 0.5),
         -0.0068433609631615745228 * (round_cna_princ6 <= 0.28661248087882996 and 
                                     round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_CSF1R <= 0.5),
          -0.041769674861334746274 * (round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5 and 
                                     round_CSF1R <= 0.5),
           0.029754318850161366738 * (round_drug_id__254 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_FGFR1 <= 0.5),
           0.069633058830960645191 * (round_cna_princ6 > 10.375198364257812),
          -0.008416674234969515056 * (round_cna_princ7 <= 14.027828216552734 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5 and 
                                     round_CSF1R <= 0.5),
           0.077766689581435774481 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_TIE2 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
          -0.066515103752374016288 * (round_drug_id__193 <= 0.5 and 
                                     round_drug_id__199 <= 0.5 and 
                                     round_drug_id__281 <= 0.5 and 
                                     round_ALK5 <= 0.5),
           -0.18941045959485877659 * (round_drug_id__199 <= 0.5 and 
                                     round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5),
           0.052688389338788568783 * (round_cna_princ6 > 0.28661248087882996 and 
                                     round_cna_princ23 > -12.51006031036377 and 
                                     round_cna_princ25 > -9.455399513244629 and 
                                     round_drug_id__35 <= 0.5),
           0.024539882222132280293 * (round_cna_princ3 > -14.671920776367188 and 
                                     round_cna_princ6 > 4.599311828613281 and 
                                     round_cna_princ18 <= 10.747390747070312 and 
                                     round_cna_princ20 > -17.991825103759766),
           0.024468900482255299955 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_FGFR1 <= 0.5),
           -0.25703697809757131143 * (round_cna_princ3 <= -14.671920776367188 and 
                                     round_drug_id__193 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5),
           -0.38270518498168631938 * (round_cna_princ3 > -14.671920776367188 and 
                                     round_drug_id__193 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5),
          -0.049658397461589058597 * (round_cna_princ18 <= 0.78175288438797 and 
                                     round_cna_princ20 <= 0.6372920274734497 and 
                                     round_cna_princ25 <= 6.179938316345215 and 
                                     round_cna_princ40 <= 19.61893081665039),
           0.033068044704310135817 * (round_cna_princ15 > -6.125102996826172 and 
                                     round_cna_princ18 <= 0.78175288438797 and 
                                     round_cna_princ25 > 6.179938316345215),
           0.069343227789812347805 * (round_347733 <= 0.5 and 
                                     round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_Num_Targets <= 2.5),
         0.00036207602496924110556 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_FGFR1 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
           0.019232188165819050218 * (round_cna_princ2 <= 0.3339686989784241 and 
                                     round_cna_princ7 <= 14.027828216552734 and 
                                     round_cna_princ27 > -10.43935775756836 and 
                                     round_cna_princ33 <= 3.91170072555542),
           0.056215406404860189526 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_FGFR1 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
           0.036589631200629614549 * (round_cna_princ18 > -1.759992003440857 and 
                                     round_drug_id__35 <= 0.5 and 
                                     round_TIE2 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
           0.096499808750121401912 * (round_cna_princ7 <= 14.027828216552734 and 
                                     round_cna_princ27 > -9.745563507080078 and 
                                     round_cna_princ33 <= 3.91170072555542 and 
                                     round_drug_id__254 <= 0.5),
          -0.062221931948676513502 * (round_cna_princ18 <= 10.747390747070312 and 
                                     round_cna_princ28 <= 5.1110520362854 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5),
           -0.10013655275709855419 * (round_cna_princ9 <= -2.331207752227783 and 
                                     round_cna_princ18 <= 0.78175288438797 and 
                                     round_cna_princ20 <= 4.372381687164307 and 
                                     round_drug_id__298 <= 0.5),
          -0.032735499981200712194 * (round_cna_princ9 <= -3.025106906890869 and 
                                     round_cna_princ18 <= 0.78175288438797 and 
                                     round_cna_princ31 > -0.7409754395484924),
            0.05723012080091735565 * (round_cna_princ7 <= 0.6541947722434998 and 
                                     round_cna_princ22 > -12.752239227294922 and 
                                     round_4633 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5),
           0.035999558532585490545 * (round_4633 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_drug_id__308 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
           0.063651296455255010343 * (round_drug_id__254 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_None > 2.5),
          -0.034717703185199975946 * (round_drug_id__193 <= 0.5 and 
                                     round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5),
           0.083490369781685364048 * (round_cna_princ4 <= 1.5352027416229248 and 
                                     round_cna_princ6 > -10.274883270263672 and 
                                     round_cna_princ13 > -5.453484535217285 and 
                                     round_cna_princ33 <= 2.605760097503662),
          -0.005318552425499407238 * (round_cna_princ13 > 2.4680399894714355 and 
                                     round_cna_princ31 > 0.39496910572052 and 
                                     round_cna_princ33 > 2.605760097503662),
          -0.041739194696625583902 * (round_cna_princ4 > 1.5352027416229248 and 
                                     round_cna_princ6 <= 1.4827982187271118 and 
                                     round_cna_princ33 <= 2.605760097503662),
           0.099224344348637694146 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_drug_id__308 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
          -0.032294467788339178971 * (round_cna_princ3 > -14.671920776367188 and 
                                     round_cna_princ23 <= -0.6465892195701599 and 
                                     round_cna_princ34 <= 3.9169163703918457),
            0.08084157418096096781 * (round_drug_id__1021 <= 0.5 and 
                                     round_FGFR1 <= 0.5 and 
                                     round_FLT3 <= 0.5 and 
                                     round_VEGFR1 <= 0.5),
          -0.023882715420441097409 * (round_drug_id__281 <= 0.5 and 
                                     round_ALK5 <= 0.5 and 
                                     round_CSF1R <= 0.5 and 
                                     round_RET <= 0.5),
           0.072043281754116339499 * (round_cna_princ7 <= 14.027828216552734 and 
                                     round_cna_princ18 > -4.536163330078125 and 
                                     round_cna_princ31 > -11.189701080322266 and 
                                     round_AFF4_mut <= 0.5),
           0.012845622674539378463 * (-4.139170169830322 < round_cna_princ5 <= 12.972505569458008 and 
                                     round_cna_princ6 > -6.260534286499023 and 
                                     round_drug_id__281 <= 0.5),
          -0.047764641312315295196 * (round_cna_princ5 <= -4.139170169830322 and 
                                     round_cna_princ20 <= 0.9346799254417419 and 
                                     round_cna_princ34 <= 0.16029495000839233 and 
                                     round_drug_id__281 <= 0.5),
           0.045621678034275088554 * (round_cna_princ2 <= 24.0517635345459 and 
                                     round_cna_princ5 > -4.139170169830322 and 
                                     round_cna_princ23 > -2.0481250286102295),
          -0.010213120866825324204 * (round_cna_princ5 > -4.139170169830322 and 
                                     round_cna_princ20 <= 5.866020679473877 and 
                                     round_cna_princ23 <= -2.0481250286102295),
          -0.067508919672434747339 * (round_cna_princ15 <= -7.65211296081543),
           0.062635975602101393034 * (round_drug_id__30 <= 0.5 and 
                                     round_drug_id__35 <= 0.5 and 
                                     round_drug_id__308 <= 0.5 and 
                                     round_VEGFR1 <= 0.5),
           0.097763837138369946711 * (round_drug_id__254 <= 0.5 and 
                                     round_drug_id__308 <= 0.5 and 
                                     round_drug_id__1023 <= 0.5 and 
                                     round_VEGFR2 <= 0.5),
          -0.094085635979659076877 * (round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5 and 
                                     round_CSF1R <= 0.5),
           -0.02861797746473799392 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK <= 0.5 and 
                                     round_CSF1R <= 0.5),
          -0.018482835706589751529 * (round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_CSF1R <= 0.5 and 
                                     round_RET <= 0.5),
          -0.054397150017849407733 * (round_cna_princ5 > -4.139170169830322 and 
                                     round_cna_princ23 <= -2.0481250286102295 and 
                                     round_cna_princ29 > 2.888566017150879),
           0.021923582208693153522 * (round_drug_id__30 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5),
          -0.040404681777428903822 * (round_drug_id__281 <= 0.5 and 
                                     round_ALK5 <= 0.5 and 
                                     round_CSF1R <= 0.5 and 
                                     round_KIT <= 0.5),
          -0.096435560147468540881 * (round_cna_princ40 <= 21.451045989990234 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1029 <= 0.5 and 
                                     round_CSF1R <= 0.5),
          -0.016274197574124498605 * (round_cna_princ18 <= -1.759992003440857 and 
                                     round_cna_princ33 <= 3.91170072555542 and 
                                     round_cna_princ40 <= 21.451045989990234 and 
                                     round_drug_id__30 <= 0.5),
           0.029176404695047838261 * (round_cna_princ13 > -6.731786251068115 and 
                                     round_cna_princ15 > -7.65211296081543 and 
                                     round_cna_princ32 > -8.759933471679688 and 
                                     round_cna_princ33 <= 2.3103909492492676),
          -0.036952142295574211261 * (round_cna_princ15 > -7.65211296081543 and 
                                     round_cna_princ31 > -4.499676704406738 and 
                                     round_cna_princ32 > -8.759933471679688 and 
                                     round_cna_princ33 > 2.3103909492492676),
          -0.030820915033321651805 * (round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_ALK5 <= 0.5 and 
                                     round_CSF1R <= 0.5),
           0.086256429915173196687 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__254 <= 0.5 and 
                                     round_drug_id__1029 <= 0.5 and 
                                     round_VEGFR <= 0.5),
           0.039423959170971696608 * (round_drug_id__30 <= 0.5 and 
                                     round_drug_id__35 <= 0.5 and 
                                     round_Num_Targets <= 7.5 and 
                                     round_VEGFR1 <= 0.5),
         -0.0051343947503073323474 * (round_drug_id__281 <= 0.5 and 
                                     round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1029 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5),
           0.018246137431808341633 * (round_cna_princ6 > 0.28661248087882996 and 
                                     round_drug_id__35 <= 0.5 and 
                                     round_drug_id__312 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5),
           -0.10093046364180793817 * (round_drug_id__298 <= 0.5 and 
                                     round_drug_id__1029 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5 and 
                                     round_CSF1R <= 0.5),
          -0.052831262615158351703 * (round_cna_princ8 > 2.6170878410339355 and 
                                     round_cna_princ17 <= 13.24669075012207 and 
                                     round_cna_princ29 <= 7.485143661499023 and 
                                     round_ARID1A_mut <= 0.5),
           0.015656592552282022623 * (round_cna_princ15 > -7.65211296081543 and 
                                     round_cna_princ17 > 0.24870412051677704 and 
                                     round_cna_princ34 > -0.14065495133399963),
          -0.036188462587825188932 * (round_cna_princ2 > -9.524709701538086 and 
                                     round_cna_princ14 > -9.997001647949219 and 
                                     round_cna_princ27 > -9.745563507080078 and 
                                     round_387266 <= 0.5),
          -0.097526958096602375803 * (round_cna_princ38 > -7.14481258392334 and 
                                     round_drug_id__281 <= 0.5 and 
                                     round_drug_id__1194 <= 0.5 and 
                                     round_CSF1R <= 0.5),
          0.0074041971628461353966 * (round_cna_princ1 <= -28.866779327392578 and 
                                     round_cna_princ10 > -22.668846130371094 and 
                                     round_cna_princ33 <= 2.3103909492492676 and 
                                     round_drug_id__298 <= 0.5),
           0.023405153937802593983 * (round_drug_id__35 <= 0.5 and 
                                     round_drug_id__1021 <= 0.5 and 
                                     round_None > 1.5 and 
                                     round_VEGFR2 <= 0.5),
           0.021571726729347477647 * (round_cna_princ10 > -22.668846130371094 and 
                                     round_cna_princ17 > 2.154477596282959 and 
                                     round_4633 <= 0.5 and 
                                     round_None > 1.5)    ])

def get_type_conversion():
    return {}
INDICATOR_COLS = []

IMPUTE_VALUES = {
    u'drug_id__1194': 0.000000,
    u'cna_princ9': -2.081475,
    u'cna_princ8': 1.462562,
    u'drug_id__298': 0.000000,
    u'cna_princ3': 3.470053,
    u'cna_princ2': 2.627788,
    u'cna_princ1': -11.863478,
    u'cna_princ7': -0.197009,
    u'cna_princ6': -0.564537,
    u'cna_princ5': -0.158841,
    u'cna_princ4': -2.239346,
    u'AFF4_mut': 0.000000,
    u'drug_id__1029': 0.000000,
    u'ALK5': 0.000000,
    u'cna_princ15': -0.781271,
    u'drug_id__6': 0.000000,
    u'MSI_FACTOR': 0.000000,
    u'drug_id__254': 0.000000,
    u'drug_id__1021': 0.000000,
    u'drug_id__312': 0.000000,
    u'drug_id__1023': 0.000000,
    u'387266': 0.000000,
    u'347733': 0.000000,
    u'VEGFR1': 0.000000,
    u'CSF1R': 0.000000,
    u'VEGFR2': 0.000000,
    u'drug_id__30': 0.000000,
    u'FLT3': 0.000000,
    u'cna_princ19': -0.810333,
    u'TP53_mut': 1.000000,
    u'2019': 0.000000,
    u'cna_princ38': -0.002936,
    u'cna_princ13': -1.524516,
    u'cna_princ34': -0.888780,
    u'drug_id__35': 0.000000,
    u'cna_princ17': 1.735889,
    u'ARID1A_mut': 0.000000,
    u'cna_princ33': 0.185144,
    u'cna_princ32': -0.321623,
    u'SRC': 0.000000,
    u'None': 7.000000,
    u'VEGFR': 0.000000,
    u'Num_Targets': 2.000000,
    u'cna_princ18': 0.306805,
    u'drug_id__281': 0.000000,
    u'RET': 0.000000,
    u'FGFR1': 0.000000,
    u'4633': 0.000000,
    u'cna_princ10': 0.974080,
    u'ALK': 0.000000,
    u'drug_id__199': 0.000000,
    u'cna_princ14': 2.927823,
    u'cna_princ31': 0.164945,
    u'1063': 0.000000,
    u'drug_id__308': 0.000000,
    u'drug_id__193': 0.000000,
    u'KIT': 0.000000,
    u'cna_princ22': -0.035392,
    u'cna_princ23': -0.892117,
    u'cna_princ20': 0.401674,
    u'cna_princ40': 0.029349,
    u'cna_princ27': -1.702273,
    u'cna_princ25': 0.083041,
    u'TIE2': 0.000000,
    u'cna_princ28': -0.197951,
    u'cna_princ29': 0.455611,}


def parse_date(x, date_format):
    """ convert date strings to numeric values. """
    try:
        # float values no longer pass isinstance(x, np.float64)
        if isinstance(x, (np.float64, float)):
            x = long_type(x)
        if '%f' in date_format and date_format.startswith('v2'):
            temp = str(x)
            if re.search('[\+-][0-9]+$', temp):
                temp = re.sub('[\+-][0-9]+$', '', temp)

            date_format = date_format[2:]
            dt = datetime.strptime(temp, date_format)
            sec = calendar.timegm(dt.timetuple())
            return sec * 1000 + dt.microsecond // 1000
        elif '%M' in date_format:
            temp = str(x)
            if re.search('[\+-][0-9]+$', temp):
                temp = re.sub('[\+-][0-9]+$', '', temp)

            return calendar.timegm(datetime.strptime(temp, date_format).timetuple())
        else:
            return datetime.strptime(str(x), date_format).toordinal()
    except:
        return float('nan')


def parse_percentage(s):
    """ remove percent sign so percentage variables can be converted to numeric """
    if isinstance(s, float):
        return s
    if isinstance(s, int):
        return float(s)
    try:
        return float(s.replace('%', ''))
    except:
        return float('nan')

def parse_nonstandard_na(s):
    """ if a column contains numbers and a unique non-numeric,
        then the non-numeric is considered to be N/A
    """
    try:
        ret = float(s)
        if np.isinf(ret):
            return float('nan')
        return ret
    except:
        return float('nan')

def parse_length(s):
    """ convert feet and inches as string to inches as numeric """
    try:
        if '"' in s and "'" in s:
            sp = s.split("'")
            return float(sp[0]) * 12 + float(sp[1].replace('"', ''))
        else:
            if "'" in s:
                return float(s.replace("'", '')) * 12
            else:
                return float(s.replace('"', ''))
    except:
        return float('nan')

def parse_currency(s):
    """ strip currency characters and commas from currency columns """
    if not isinstance(s, text_type):
        return float('nan')
    s = re.sub(u'[\$\u20AC\u00A3\uFFE1\u00A5\uFFE5]|(EUR)', '', s)
    s = s.replace(',', '')
    try:
        return float(s)
    except:
        return float('nan')


def parse_currency_replace_cents_period(val, currency_symbol):
    try:
        if np.isnan(val):
            return val
    except TypeError:
        pass
    if not isinstance(val, string_types):
        raise ValueError('Found wrong value for currency: {}'.format(val))
    try:
        val = val.replace(currency_symbol, "", 1)
        val = val.replace(" ", "")
        val = val.replace(",", "")
        val = float(val)
    except ValueError:
        val = float('nan')
    return val


def parse_currency_replace_cents_comma(val, currency_symbol):
    try:
        if np.isnan(val):
            return val
    except TypeError:
        pass
    if not isinstance(val, string_types):
        raise ValueError('Found wrong value for currency: {}'.format(val))
    try:
        val = val.replace(currency_symbol, "", 1)
        val = val.replace(" ", "")
        val = val.replace(".", "")
        val = val.replace(",", ".")
        val = float(val)
    except ValueError:
        val = float('nan')
    return val


def parse_currency_replace_no_cents(val, currency_symbol):
    try:
        if np.isnan(val):
            return val
    except TypeError:
        pass
    if not isinstance(val, string_types):
        raise ValueError('Found wrong value for currency: {}'.format(val))
    try:
        val = val.replace(currency_symbol, "", 1)
        val = val.replace(" ", "")
        val = val.replace(",", "")
        val = val.replace(".", "")
        val = float(val)
    except ValueError:
        val = float('nan')
    return val

def parse_numeric_types(ds):
    """ convert strings with numeric types (date, currency, etc.)
        to actual numeric values """
    TYPE_CONVERSION = get_type_conversion()
    for col in ds.columns:
        if col in TYPE_CONVERSION:
            convert_func = TYPE_CONVERSION[col]['convert_func']
            convert_args = TYPE_CONVERSION[col]['convert_args']
            ds[col] = ds[col].apply(convert_func, args=convert_args)
    return ds

def sanitize_name(name):
    safe = name.strip().replace("-", "_").replace("$", "_").replace(".", "_")
    safe = safe.replace("{", "_").replace("}", "_")
    safe = safe.replace('"', '_')
    return safe

def rename_columns(ds):
    new_names = {}
    existing_names = set()
    disambiguation = {}
    blank_index = 0
    for old_col in ds.columns:
        col = sanitize_name(old_col)
        if col == '':
            col = 'Unnamed: %d' % blank_index
            blank_index += 1
        if col in existing_names:
            suffix = '_%d' % disambiguation.setdefault(col, 1)
            disambiguation[col] += 1
            col = col + suffix
        existing_names.add(col)
        new_names[old_col] = col
    ds.rename(columns=new_names, inplace=True)
    return ds

def add_missing_indicators(ds):
    for col in INDICATOR_COLS:
        ds[col + '-mi'] = ds[col].isnull().astype(int)
    return ds

def impute_values(ds):
    for col in ds:
        if col in IMPUTE_VALUES:
            ds.loc[ds[col].isnull(), col] = IMPUTE_VALUES[col]
    return ds

BIG_LEVELS = {
}


SMALL_NULLS = {
}


VAR_TYPES = {
    u'drug_id__1194': 'N',
    u'cna_princ9': 'N',
    u'cna_princ8': 'N',
    u'drug_id__298': 'N',
    u'cna_princ3': 'N',
    u'cna_princ2': 'N',
    u'cna_princ1': 'N',
    u'cna_princ7': 'N',
    u'cna_princ6': 'N',
    u'cna_princ5': 'N',
    u'cna_princ4': 'N',
    u'AFF4_mut': 'N',
    u'drug_id__1029': 'N',
    u'ALK5': 'N',
    u'cna_princ15': 'N',
    u'drug_id__6': 'N',
    u'MSI_FACTOR': 'N',
    u'drug_id__254': 'N',
    u'drug_id__1021': 'N',
    u'drug_id__312': 'N',
    u'drug_id__1023': 'N',
    u'387266': 'N',
    u'347733': 'N',
    u'VEGFR1': 'N',
    u'CSF1R': 'N',
    u'VEGFR2': 'N',
    u'drug_id__30': 'N',
    u'FLT3': 'N',
    u'cna_princ19': 'N',
    u'TP53_mut': 'N',
    u'2019': 'N',
    u'cna_princ38': 'N',
    u'cna_princ13': 'N',
    u'cna_princ34': 'N',
    u'drug_id__35': 'N',
    u'cna_princ31': 'N',
    u'ARID1A_mut': 'N',
    u'cna_princ33': 'N',
    u'cna_princ32': 'N',
    u'SRC': 'N',
    u'None': 'N',
    u'VEGFR': 'N',
    u'Num_Targets': 'N',
    u'cna_princ18': 'N',
    u'drug_id__281': 'N',
    u'RET': 'N',
    u'FGFR1': 'N',
    u'4633': 'N',
    u'cna_princ10': 'N',
    u'ALK': 'N',
    u'drug_id__199': 'N',
    u'cna_princ14': 'N',
    u'cna_princ17': 'N',
    u'1063': 'N',
    u'drug_id__308': 'N',
    u'drug_id__193': 'N',
    u'KIT': 'N',
    u'cna_princ22': 'N',
    u'cna_princ23': 'N',
    u'cna_princ20': 'N',
    u'cna_princ40': 'N',
    u'cna_princ27': 'N',
    u'cna_princ25': 'N',
    u'TIE2': 'N',
    u'cna_princ28': 'N',
    u'cna_princ29': 'N',
}


def combine_small_levels(ds):
    for col in ds:
        if BIG_LEVELS.get(col, None) is not None:
            mask = np.logical_and(~ds[col].isin(BIG_LEVELS[col]), ds[col].notnull())
            if np.any(mask):
                ds.loc[mask, col] = 'small_count'
        if SMALL_NULLS.get(col):
            mask = ds[col].isnull()
            if np.any(mask):
                ds.loc[mask, col] = 'small_count'
        if VAR_TYPES.get(col) == 'C' or VAR_TYPES.get(col) == 'T':
            mask = ds[col].isnull()
            if np.any(mask):
                if ds[col].dtype == float:
                    ds[col] = ds[col].astype(object)
                ds.loc[mask, col] = 'nan'
    return ds

# N/A strings in addition to the ones used by Pandas read_csv()
NA_VALUES = ['null', 'na', 'n/a', '#N/A', 'N/A', '?', '.', '', 'Inf', 'INF', 'inf', '-inf', '-Inf', '-INF', ' ', 'None', 'NaN', '-nan', 'NULL', 'NA', '-1.#IND', '1.#IND', '-1.#QNAN', '1.#QNAN', '#NA', '#N/A N/A', '-NaN', 'nan']

# True/False strings in addition to the ones used by Pandas read_csv()
TRUE_VALUES = ['TRUE', 'True', 'true']
FALSE_VALUES = ['FALSE', 'False', 'false']

DEFAULT_ENCODING = 'utf8'

REQUIRED_COLUMNS = [u"drug_id__1194",u"cna_princ9",u"ALK5",u"drug_id__298",u"cna_princ3",u"cna_princ2",u"cna_princ1",u"cna_princ7",u"cna_princ6",u"RET",u"cna_princ4",u"AFF4_mut",u"cna_princ8",u"cna_princ10",u"drug_id__6",u"MSI_FACTOR",u"drug_id__254",u"drug_id__1021",u"drug_id__312",u"drug_id__1023",u"387266",u"347733",u"drug_id__1029",u"CSF1R",u"cna_princ20",u"FLT3",u"cna_princ19",u"cna_princ18",u"2019",u"cna_princ38",u"cna_princ13",u"TP53_mut",u"drug_id__35",u"cna_princ31",u"ARID1A_mut",u"drug_id__30",u"cna_princ27",u"SRC",u"None",u"cna_princ33",u"VEGFR",u"Num_Targets",u"drug_id__281",u"cna_princ5",u"cna_princ34",u"cna_princ32",u"FGFR1",u"4633",u"cna_princ29",u"ALK",u"drug_id__199",u"cna_princ14",u"cna_princ17",u"1063",u"drug_id__308",u"drug_id__193",u"KIT",u"cna_princ22",u"cna_princ23",u"VEGFR2",u"cna_princ40",u"cna_princ15",u"cna_princ25",u"TIE2",u"cna_princ28",u"VEGFR1"]


def validate_columns(column_list):
    if set(REQUIRED_COLUMNS) <= set(column_list):
        return True
    else :
        raise ValueError("Required columns missing: %s" %
                         (set(REQUIRED_COLUMNS) - set(column_list)))

def convert_bool(ds):
    TYPE_CONVERSION = get_type_conversion()
    for col in ds.columns:
        if VAR_TYPES.get(col) == 'C' and ds[col].dtype in (int, float):
            mask = ds[col].notnull()
            ds[col] = ds[col].astype(object)
            ds.loc[mask, col] = ds.loc[mask, col].astype(text_type)
        elif VAR_TYPES.get(col) == 'N' and ds[col].dtype == bool:
            ds[col] = ds[col].astype(float)
        elif ds[col].dtype == bool:
            ds[col] = ds[col].astype(text_type)
        elif ds[col].dtype == object:
            if VAR_TYPES.get(col) == 'N' and col not in TYPE_CONVERSION:
                mask = ds[col].apply(lambda x: x in TRUE_VALUES)
                if np.any(mask):
                    ds.loc[mask, col] = 1
                mask = ds[col].apply(lambda x: x in FALSE_VALUES)
                if np.any(mask):
                    ds.loc[mask, col] = 0
                ds[col] = ds[col].astype(float)
            elif TYPE_CONVERSION.get(col) is None:
                mask = ds[col].notnull()
                ds.loc[mask, col] = ds.loc[mask, col].astype(text_type)
    return ds

def get_dtypes():
    return {a: object for a, b in VAR_TYPES.items() if b == 'C'}

def predict_dataframe(ds):
    return ds.apply(predict, axis=1)

def run_dataframe(ds):
    ds = rename_columns(ds)
    ds = convert_bool(ds)
    validate_columns(ds.columns)
    ds = parse_numeric_types(ds)
    ds = add_missing_indicators(ds)
    ds = impute_values(ds)
    ds = combine_small_levels(ds)
    prediction = predict_dataframe(ds)
    return prediction


def run(dataset_path, output_path, encoding=None):
    if encoding is None:
        encoding = DEFAULT_ENCODING

    ds = pd.read_csv(dataset_path, na_values=NA_VALUES, low_memory=False,
                     dtype=get_dtypes(), encoding=encoding)

    prediction = run_dataframe(ds)
    prediction_file = output_path
    prediction.name = 'Prediction'
    prediction.to_csv(prediction_file, header=True, index_label='Index')


def _construct_parser():
    import argparse

    parser = argparse.ArgumentParser(description='Make offline predictions with DataRobot Prime')

    parser.add_argument(
        '--encoding',
        type=str,
        help=('the encoding of the dataset you are going to make predictions with. '
              'DataRobot Prime defaults to UTF-8 if not otherwise specified. See the '
              '"Codecs" column of the Python-supported standards chart '
              '(https://docs.python.org/2/library/codecs.html#standard-encodings) '
              'for possible alternative entries.'),
        metavar='<encoding>'
    )
    parser.add_argument(
        'input_path',
        type=str,
        help=('a .csv file (your dataset); columns must correspond to the '
              'feature set used to generate the DataRobot Prime model.'),
        metavar='<data_file>'
    )
    parser.add_argument(
        'output_path',
        type=str,
        help='the filename where DataRobot writes the results.',
        metavar='<output_file>'
    )

    return parser


def _parse_command(args):
    parser = _construct_parser()
    parsed_args = parser.parse_args(args[1:])

    if parsed_args.encoding is None:
        sys.stderr.write('Warning: For input data encodings other than the standard utf-8, '
                         'see documentation at http://10.132.37.7/docs/users-guide/more-info/tabs/prime-examples.html')
        parsed_args.encoding = DEFAULT_ENCODING

    return parsed_args


if __name__ == '__main__':
    args = _parse_command(sys.argv)
    run(args.input_path, args.output_path, encoding=args.encoding)
