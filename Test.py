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

RANK_FNAME = "Cas9_score_results/RANK_final_DeepCas9_Final_"

SORTED_EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_w_cas9_scores_"

PAM_SEQ = ['NGG']
ADD_SEQ1_LEN = 4
SPACER_LEN = 20
ADD_SEQ2_LEN = 3
CLVG_AFTER_PAM = 3

FILE_NAME_LIST = ["X", "Y"]

INITIAL_SEQ = [PAM_SEQ, ADD_SEQ1_LEN, SPACER_LEN, ADD_SEQ2_LEN, CLVG_AFTER_PAM, WORK_DIR]

############### CAS OFF FINDER ################
CAS_OFF_EXCEL = "Chlorocebus_sabaeus_excel_20200424/Chlorocebus_sabaeus_filter_out_05_65_20200506/chlorochebus_sabaeus_w_cas9_scores_"
CAS_OFF_TXT = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_cas_NGG_off"
MAX_MISMATCH = 3
SUB_SET_NUM = 100

INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, SUB_SET_NUM]
############### end setting env ################

def main1():
    util = Util.Utils()
    logic = Logic.Logics()
    logic_pre = LogicPrep.LogicsPrep()
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


    # df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, "Y")
    # data_obj = df_obj.to_dict()
    # query_seq_dict = data_obj["order sgRNA Target sequence"]
    # pam_dict = data_obj["order PAM"]
    #
    # for query_idx in range(len(query_seq_dict)):
    #     cas_input_set.add(query_seq_dict[query_idx] + pam_dict[query_idx])

    print("total dict len : " + str(total_len))
    print(" ")
    print(" ")

    util.make_txt_cas_off_finder_input(WORK_DIR + CAS_OFF_TXT, cas_input_set, INITIAL_CAS_OFF)






def main2():
    num_key_dict = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
    for key, val in num_key_dict.items():
        print(str(key) + ":" + str(val))
    pass

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# main_merge()
make_cas_off_finder_input()
# main1()
# main2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))