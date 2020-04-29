from time import clock
import re
import numpy as np


import Util
import Logic
import LogicPre
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/SongMyunJae_YuGooSang/20200417/WORK_DIR/"

SEQ_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome."

CDS_EACH_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chromosome."
CDS_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chr"

EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_"

RANK_FNAME = "Cas9_score_results/RANK_final_DeepCas9_Final_"

SORTED_EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_w_cas9_scores_"

PAM_SEQ = ['NGG']
ADD_SEQ1_LEN = 4
SPACER_LEN = 20
ADD_SEQ2_LEN = 3
CLVG_AFTER_PAM = 3

FILE_NAME_LIST = ["X", "Y"]


INITIAL_SEQ = [PAM_SEQ, ADD_SEQ1_LEN, SPACER_LEN, ADD_SEQ2_LEN, CLVG_AFTER_PAM, WORK_DIR]
############### end setting env  ################

def main1():
    util = Util.Utils()
    logic = Logic.Logics()
    logic_pre = LogicPre.LogicsPre()
    valid = Valid.Validations()


def main_merge():
    util = Util.Utils()
    valid = Valid.Validations()

    for i in range(1,30):
        FILE_NAME_LIST.append(str(i))

    for f_num in FILE_NAME_LIST:
        tmp_arr = util.read_cas9_scores_arr(WORK_DIR + RANK_FNAME, f_num)
        df_obj = util.read_excel_2_dataframe(WORK_DIR + EXCEL_FNAME, f_num)
        if valid.is_same_len(tmp_arr, df_obj, f_num):
            continue

        tmp_dict = util.merge_data(tmp_arr, df_obj)

        print(f_num + " is DONE")
        util.make_excel_srt_by_deepcas9_scr(WORK_DIR + SORTED_EXCEL_FNAME, f_num, tmp_dict)

def main2():
    num_key_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    for key, val in num_key_dict.items():
        print(str(key) + ":" + str(val))
    pass

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
main_merge()
# main1()
# main2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))