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
    round_Antimetabolite = np.float32(row[u'Antimetabolite'])
    round_322 = np.float32(row[u'322'])
    round_ASPM_mut = np.float32(row[u'ASPM_mut'])
    round_790 = np.float32(row[u'790'])
    round_cna_princ12 = np.float32(row[u'cna_princ12'])
    round_cna_princ49 = np.float32(row[u'cna_princ49'])
    round_cna_princ48 = np.float32(row[u'cna_princ48'])
    round_cna_princ47 = np.float32(row[u'cna_princ47'])
    round_cna_princ46 = np.float32(row[u'cna_princ46'])
    round_DNAcrosslinker = np.float32(row[u'DNAcrosslinker'])
    round_ASXL1_mut = np.float32(row[u'ASXL1_mut'])
    round_cna_princ40 = np.float32(row[u'cna_princ40'])
    round_cna_princ37 = np.float32(row[u'cna_princ37'])
    round_drug_id__1494 = np.float32(row[u'drug_id__1494'])
    round_cna_princ18 = np.float32(row[u'cna_princ18'])
    round_cna_princ33 = np.float32(row[u'cna_princ33'])
    round_cna_princ30 = np.float32(row[u'cna_princ30'])
    round_cna_princ43 = np.float32(row[u'cna_princ43'])
    round_cna_princ16 = np.float32(row[u'cna_princ16'])
    round_cna_princ17 = np.float32(row[u'cna_princ17'])
    round_cna_princ10 = np.float32(row[u'cna_princ10'])
    round_cna_princ38 = np.float32(row[u'cna_princ38'])
    round_cna_princ13 = np.float32(row[u'cna_princ13'])
    round_drug_id__1005 = np.float32(row[u'drug_id__1005'])
    round_10075 = np.float32(row[u'10075'])
    round_ARID4B_mut = np.float32(row[u'ARID4B_mut'])
    round_drug_id__190 = np.float32(row[u'drug_id__190'])
    round_drug_id__135 = np.float32(row[u'drug_id__135'])
    round_drug_id__134 = np.float32(row[u'drug_id__134'])
    round_drug_id__136 = np.float32(row[u'drug_id__136'])
    round_dsDNAbreakinduction = np.float32(row[u'dsDNAbreakinduction'])
    round_drug_id__133 = np.float32(row[u'drug_id__133'])
    round_cna_princ19 = np.float32(row[u'cna_princ19'])
    round_cna_princ25 = np.float32(row[u'cna_princ25'])
    round_cna_princ27 = np.float32(row[u'cna_princ27'])
    round_cna_princ26 = np.float32(row[u'cna_princ26'])
    round_cna_princ21 = np.float32(row[u'cna_princ21'])
    round_cna_princ20 = np.float32(row[u'cna_princ20'])
    round_drug_id__1375 = np.float32(row[u'drug_id__1375'])
    round_cna_princ22 = np.float32(row[u'cna_princ22'])
    round_cna_princ28 = np.float32(row[u'cna_princ28'])
    round_51168 = np.float32(row[u'51168'])
    round_cna_princ9 = np.float32(row[u'cna_princ9'])
    round_cna_princ7 = np.float32(row[u'cna_princ7'])
    round_cna_princ4 = np.float32(row[u'cna_princ4'])
    round_cna_princ5 = np.float32(row[u'cna_princ5'])
    round_cna_princ3 = np.float32(row[u'cna_princ3'])
    round_TP53_mut = np.float32(row[u'TP53_mut'])
    round_cna_princ1 = np.float32(row[u'cna_princ1'])
    return sum([
         0.0034698,
           -0.10396085248440986581 * (round_drug_id__135),
             2.5598146299306492324 * (round_drug_id__1375),
             -1.827783443069198599 * (round_drug_id__1494),
            0.22754105240382252462 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           -0.21729765432472944453 * (round_cna_princ9 <= 6.698186874389648 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
          -0.076999144084847662151 * (round_cna_princ3 <= 24.60089111328125 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
           0.036145415648419727817 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           -0.27570165320657191721 * (round_cna_princ9 <= 25.0417537689209 and 
                                     round_drug_id__190 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
         -0.0091857237112148680536 * (round_cna_princ10 > -2.2314062118530273 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
            0.56623720988132830723 * (round_790 <= 0.5 and 
                                     round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
         -0.0042434029505440153118 * (round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
           -0.23577007065159735766 * (round_cna_princ4 <= 26.178735733032227 and 
                                     round_cna_princ13 <= -0.7604091167449951 and 
                                     round_cna_princ33 <= 3.0836405754089355 and 
                                     round_drug_id__1005 <= 0.5),
          -0.039457107328061513662 * (round_cna_princ13 <= -4.665262222290039 and 
                                     round_cna_princ17 > -10.778310775756836 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
           -0.29284951553160443183 * (round_cna_princ9 <= 25.0417537689209 and 
                                     round_cna_princ48 <= 4.09360408782959 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
          -0.016516609944656524683 * (round_cna_princ26 <= 6.804866313934326 and 
                                     round_cna_princ37 > 5.084649085998535 and 
                                     round_drug_id__1375 <= 0.5),
           -0.33260864199378264994 * (round_cna_princ10 > -2.2314062118530273 and 
                                     round_cna_princ13 <= 21.811073303222656 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_ASXL1_mut <= 0.5),
            0.67349622916665219385 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5 and 
                                     round_Antimetabolite <= 0.5),
            0.10023249532720066024 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           -0.11137879381779396126 * (round_cna_princ7 <= 14.678129196166992 and 
                                     round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_ASXL1_mut <= 0.5),
           0.061529091722015608057 * (round_cna_princ13 > -0.7604091167449951 and 
                                     round_cna_princ21 > -4.157510757446289 and 
                                     round_cna_princ27 > 1.826384425163269),
          -0.094212291135550005339 * (round_cna_princ26 <= 4.672094345092773 and 
                                     round_cna_princ30 > -4.092028617858887 and 
                                     round_cna_princ43 <= 2.1591529846191406 and 
                                     round_ASXL1_mut <= 0.5),
           0.015444782866292257328 * (round_cna_princ10 <= 12.115333557128906 and 
                                     round_cna_princ43 > 2.1591529846191406 and 
                                     round_cna_princ49 > -3.859037399291992 and 
                                     round_ASXL1_mut <= 0.5),
            0.13848340117300997343 * (round_cna_princ4 > -4.572201728820801 and 
                                     round_cna_princ28 > 4.374505043029785 and 
                                     round_cna_princ46 <= 1.9285870790481567),
           -0.16715574544796621415 * (round_cna_princ1 > -26.582984924316406 and 
                                     round_cna_princ4 <= -4.572201728820801 and 
                                     round_drug_id__1375 <= 0.5),
          -0.062375393632863804227 * (round_cna_princ9 <= 6.698186874389648 and 
                                     round_cna_princ47 <= -1.2316865921020508 and 
                                     round_cna_princ48 <= 4.09360408782959 and 
                                     round_drug_id__1375 <= 0.5),
            0.18713613067916143318 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
          -0.022228401214528670932 * (round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_DNAcrosslinker <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
            0.14511708949284385373 * (round_cna_princ38 > -6.957921981811523 and 
                                     round_51168 <= 0.5 and 
                                     round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5),
          -0.032758632992544500784 * (round_cna_princ48 <= 4.09360408782959 and 
                                     round_drug_id__190 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
           -0.16108532129329578009 * (round_cna_princ20 > -15.553276062011719 and 
                                     round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
         -0.0031926831268411478074 * (round_cna_princ3 > -10.972679138183594 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
          0.0098135461976269471734 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
          -0.063321263563560478693 * (round_cna_princ7 <= 16.94678497314453 and 
                                     round_cna_princ25 > -5.701957702636719 and 
                                     round_cna_princ43 <= 2.1591529846191406 and 
                                     round_cna_princ48 <= 3.0994534492492676),
          -0.073018059085696973565 * (round_cna_princ3 > -10.972679138183594 and 
                                     round_cna_princ9 <= 25.0417537689209 and 
                                     round_cna_princ37 > 3.0032124519348145 and 
                                     round_drug_id__1005 <= 0.5),
            0.15975393055993217239 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           -0.18332959478759208727 * (round_cna_princ5 <= 13.567386627197266 and 
                                     round_cna_princ43 <= 2.1591529846191406 and 
                                     round_cna_princ47 <= -0.5378657579421997 and 
                                     round_ASPM_mut <= 0.5),
           -0.20261389928089790224 * (round_cna_princ16 > -2.7181954383850098 and 
                                     round_cna_princ20 <= 7.9051313400268555 and 
                                     round_cna_princ43 > -5.000893592834473 and 
                                     round_drug_id__1005 <= 0.5),
        -2.3246732970774360569E-05 * (round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
           0.037600666121716544432 * (round_cna_princ13 > -0.7604091167449951 and 
                                     round_cna_princ37 > -1.5832033157348633 and 
                                     round_cna_princ48 > -8.795920372009277),
           0.056485378183827615606 * (round_cna_princ10 <= -2.2314062118530273 and 
                                     round_cna_princ12 <= 10.9129056930542 and 
                                     round_drug_id__133 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
          -0.019134766609728303538 * (round_cna_princ3 > -10.972679138183594 and 
                                     round_cna_princ5 > -21.135290145874023 and 
                                     round_cna_princ40 > -4.627498626708984 and 
                                     round_10075 <= 0.5),
           0.079763434605173433312 * (round_cna_princ37 <= 3.0032124519348145 and 
                                     round_cna_princ46 > -0.20845897495746613),
           0.074456246414085769381 * (round_cna_princ18 > 0.26090186834335327 and 
                                     round_cna_princ40 <= 6.644059181213379 and 
                                     round_ARID4B_mut <= 0.5 and 
                                     round_ASXL1_mut <= 0.5),
           0.030641572315460295639 * (round_cna_princ10 <= 12.115333557128906 and 
                                     2.1591529846191406 < round_cna_princ43 <= 6.772865295410156 and 
                                     round_ASPM_mut <= 0.5),
          -0.086158488897822502195 * (round_cna_princ3 > -10.972679138183594 and 
                                     round_cna_princ7 <= 14.678129196166992 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_ASPM_mut <= 0.5),
            0.45924461902344165809 * (round_322 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           0.011094920835491024758 * (round_cna_princ19 > -9.839859962463379 and 
                                     round_cna_princ37 <= 3.0032124519348145 and 
                                     round_10075 <= 0.5),
            0.13440148383907057905 * (round_drug_id__133 <= 0.5 and 
                                     round_drug_id__135 <= 0.5 and 
                                     round_drug_id__136 <= 0.5 and 
                                     round_drug_id__1494 <= 0.5),
           -0.50188414017717875026 * (round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5 and 
                                     round_dsDNAbreakinduction <= 0.5),
           -0.12607139325618582637 * (round_cna_princ25 > -14.158241271972656 and 
                                     round_drug_id__134 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5),
          0.0039105402674197188573 * (round_cna_princ37 <= 3.0032124519348145 and 
                                     round_cna_princ46 > -0.20845897495746613 and 
                                     round_drug_id__134 <= 0.5),
            0.13001850422503108407 * (round_cna_princ1 <= 85.21183776855469 and 
                                     round_cna_princ4 > 1.04206120967865 and 
                                     round_cna_princ38 > -4.000243186950684 and 
                                     round_TP53_mut > 0.5),
          -0.063620242676072108479 * (round_cna_princ22 <= 23.697113037109375 and 
                                     round_drug_id__190 <= 0.5 and 
                                     round_drug_id__1005 <= 0.5 and 
                                     round_drug_id__1375 <= 0.5)    ])

def get_type_conversion():
    return {}
INDICATOR_COLS = []

IMPUTE_VALUES = {
    u'drug_id__1375': 0.000000,
    u'cna_princ9': -3.083993,
    u'dsDNAbreakinduction': 0.000000,
    u'cna_princ3': -2.046762,
    u'cna_princ1': -11.040401,
    u'cna_princ7': -1.660806,
    u'cna_princ5': 1.746064,
    u'cna_princ4': -0.758737,
    u'TP53_mut': 1.000000,
    u'DNAcrosslinker': 0.000000,
    u'790': 0.000000,
    u'cna_princ46': 0.209352,
    u'drug_id__1005': 0.000000,
    u'322': 0.000000,
    u'ASPM_mut': 0.000000,
    u'cna_princ19': -0.996765,
    u'cna_princ47': -0.448365,
    u'cna_princ38': 1.486623,
    u'cna_princ13': -0.622516,
    u'cna_princ40': -0.458685,
    u'cna_princ37': 1.166323,
    u'cna_princ10': 1.478655,
    u'cna_princ17': -1.035451,
    u'cna_princ30': -0.579028,
    u'cna_princ33': -0.617424,
    u'10075': 0.000000,
    u'cna_princ43': -0.644742,
    u'ARID4B_mut': 0.000000,
    u'Antimetabolite': 0.000000,
    u'cna_princ12': 2.567354,
    u'drug_id__133': 0.000000,
    u'drug_id__134': 0.000000,
    u'drug_id__135': 0.000000,
    u'drug_id__136': 0.000000,
    u'cna_princ18': -1.414706,
    u'51168': 0.000000,
    u'drug_id__190': 0.000000,
    u'cna_princ16': 0.239906,
    u'cna_princ22': -1.051865,
    u'drug_id__1494': 0.000000,
    u'cna_princ20': -1.328037,
    u'cna_princ21': -2.445388,
    u'cna_princ26': -0.019441,
    u'cna_princ27': 0.687658,
    u'cna_princ25': -0.793549,
    u'cna_princ28': 0.636345,
    u'ASXL1_mut': 0.000000,
    u'cna_princ48': 0.044565,
    u'cna_princ49': -0.472128,}


def bag_of_words(text):
    """ set of whole words  in a block of text """
    if type(text) == float:
        return set()

    return set(word.lower() for word in
               re.findall(r'\w+', text, re.UNICODE | re.IGNORECASE))


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
    u'drug_id__1375': 'N',
    u'cna_princ9': 'N',
    u'dsDNAbreakinduction': 'N',
    u'cna_princ3': 'N',
    u'cna_princ1': 'N',
    u'cna_princ7': 'N',
    u'cna_princ5': 'N',
    u'cna_princ4': 'N',
    u'TP53_mut': 'N',
    u'DNAcrosslinker': 'N',
    u'790': 'N',
    u'cna_princ46': 'N',
    u'drug_id__1005': 'N',
    u'322': 'N',
    u'ASPM_mut': 'N',
    u'cna_princ19': 'N',
    u'cna_princ18': 'N',
    u'cna_princ38': 'N',
    u'cna_princ13': 'N',
    u'cna_princ40': 'N',
    u'cna_princ37': 'N',
    u'cna_princ10': 'N',
    u'cna_princ17': 'N',
    u'cna_princ30': 'N',
    u'cna_princ33': 'N',
    u'10075': 'N',
    u'cna_princ43': 'N',
    u'ARID4B_mut': 'N',
    u'Antimetabolite': 'N',
    u'cna_princ12': 'N',
    u'drug_id__133': 'N',
    u'drug_id__134': 'N',
    u'drug_id__135': 'N',
    u'drug_id__136': 'N',
    u'cna_princ47': 'N',
    u'51168': 'N',
    u'drug_id__190': 'N',
    u'cna_princ16': 'N',
    u'cna_princ22': 'N',
    u'drug_id__1494': 'N',
    u'cna_princ20': 'N',
    u'cna_princ21': 'N',
    u'cna_princ26': 'N',
    u'cna_princ27': 'N',
    u'cna_princ25': 'N',
    u'cna_princ28': 'N',
    u'ASXL1_mut': 'N',
    u'cna_princ48': 'N',
    u'cna_princ49': 'N',
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

REQUIRED_COLUMNS = [u"drug_id__1375",u"cna_princ9",u"dsDNAbreakinduction",u"cna_princ3",u"cna_princ1",u"cna_princ7",u"cna_princ5",u"cna_princ4",u"TP53_mut",u"DNAcrosslinker",u"790",u"cna_princ46",u"drug_id__1005",u"322",u"ASPM_mut",u"cna_princ19",u"cna_princ47",u"cna_princ38",u"cna_princ13",u"cna_princ12",u"cna_princ37",u"cna_princ10",u"cna_princ17",u"cna_princ30",u"cna_princ33",u"10075",u"cna_princ40",u"cna_princ43",u"cna_princ18",u"ARID4B_mut",u"Antimetabolite",u"drug_id__133",u"drug_id__134",u"drug_id__135",u"drug_id__136",u"51168",u"drug_id__190",u"cna_princ16",u"cna_princ22",u"drug_id__1494",u"cna_princ20",u"cna_princ21",u"cna_princ26",u"cna_princ27",u"cna_princ25",u"cna_princ28",u"ASXL1_mut",u"cna_princ48",u"cna_princ49"]


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
