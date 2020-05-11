from time import clock
import re
import numpy as np


import Util
import Logic
import LogicPrep
import Valid
############### start to set env ################
WORK_DIR = "D:/000_WORK/SongMyunJae_YuGooSang/20200417/WORK_DIR/"

SEQ_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome."

CDS_EACH_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chromosome."
CDS_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chr"

EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_"

PAM_SEQ = ['NGG']
ADD_SEQ1_LEN = 4
SPACER_LEN = 20
ADD_SEQ2_LEN = 3
CLVG_AFTER_PAM = 3

FILE_NAME_LIST = ["X", "Y"]

INITIAL_SEQ = [PAM_SEQ, ADD_SEQ1_LEN, SPACER_LEN, ADD_SEQ2_LEN, CLVG_AFTER_PAM, WORK_DIR]

############### MERGE DEEP_CAS_9 ##############
RANK_FNAME = "Cas9_score_results/RANK_final_DeepCas9_Final_"
SORTED_EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_w_cas9_scores_"

############### CAS OFF FINDER ################
CAS_OFF_EXCEL = "Chlorocebus_sabaeus_excel_20200424/Chlorocebus_sabaeus_filter_out_05_65_20200506/chlorochebus_sabaeus_w_cas9_scores_"
CAS_OFF_TXT = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_cas_NGG_off_"
MAX_MISMATCH = 3
SUB_SET_NUM = 100
REF_SRV_PATH = "chlorocebus_sabaeus_chr"

INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, SUB_SET_NUM, WORK_DIR + CAS_OFF_TXT, REF_SRV_PATH]

########### MERGE OFF_TAGET RESULT #############
CNT_RESULT_PATH = "CAS_OFF_FINDER/Output/CountResult/chlorochebus_sabaeus_cas_NGG_off"
CNT_RESULT_EXT = "_result_count.txt"
FILE_CNT = 34

INITIAL_MRGE_OFF_TRGT = [WORK_DIR + CNT_RESULT_PATH, CNT_RESULT_EXT, FILE_CNT]
############### end setting env ################

def merge_off_target():
    util = Util.Utils()
    logic = Logic.Logics()
    logic_pre = LogicPrep.LogicsPrep()
    valid = Valid.Validations()

    off_trgt_dict = util.read_txt_to_dict(WORK_DIR + CNT_RESULT_PATH, CNT_RESULT_EXT, FILE_CNT)

    # for i in range(1, 30):
    #     FILE_NAME_LIST.append(str(i))

    for f_num in FILE_NAME_LIST:
        print(str(f_num))
        df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, f_num)

        util.make_excel_off_target_data(off_trgt_dict, df_obj, WORK_DIR + CAS_OFF_EXCEL + "off_trgt_", f_num)

def merge_deep_cas_9():
    util = Util.Utils()
    valid = Valid.Validations()

    for i in range(1,30):
        FILE_NAME_LIST.append(str(i))

    for f_num in FILE_NAME_LIST:
        tmp_arr = util.read_cas9_scores_arr(WORK_DIR + RANK_FNAME, f_num)
        df_obj = util.read_excel_2_dataframe(WORK_DIR + EXCEL_FNAME, f_num)
        if valid.is_same_len(tmp_arr, df_obj, f_num):
            continue

        tmp_dict = util.merge_deep_cas_9_data(tmp_arr, df_obj)

        print(f_num + " is DONE")
        util.make_excel_srt_by_deepcas9_scr(WORK_DIR + SORTED_EXCEL_FNAME, f_num, tmp_dict)

def make_cas_off_finder_input():
    util = Util.Utils()
    cas_input_set = set()
    for i in range(1,30):
        FILE_NAME_LIST.append(str(i))

    total_len = 0
    for f_num in FILE_NAME_LIST:
        df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, f_num)

        data_obj = df_obj.to_dict()
        query_seq_dict = data_obj["order sgRNA Target sequence"]
        pam_dict = data_obj["order PAM"]
        query_seq_dict_len = len(query_seq_dict)
        total_len += query_seq_dict_len
        print("file [" + str(f_num) + "] len : " + str(query_seq_dict_len))

        for query_idx in range(query_seq_dict_len):
            cas_input_set.add(query_seq_dict[query_idx] + pam_dict[query_idx])

    util.make_txt_cas_off_finder_input(cas_input_set, INITIAL_CAS_OFF)






def main2():
    num_key_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    for key, val in num_key_dict.items():
        print(str(key) + ":" + str(val))
    pass

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
merge_off_target()
# merge_deep_cas_9()
# make_cas_off_finder_input()
# main2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))