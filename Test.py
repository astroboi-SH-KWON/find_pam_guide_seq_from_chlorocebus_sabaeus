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
# FILE_NAME_LIST = ["Y"]

INITIAL_SEQ = [PAM_SEQ, ADD_SEQ1_LEN, SPACER_LEN, ADD_SEQ2_LEN, CLVG_AFTER_PAM, WORK_DIR]

############### MERGE DEEP_CAS_9 ##############
RANK_FNAME = "Cas9_score_results/RANK_final_DeepCas9_Final_"
SORTED_EXCEL_FNAME = "Chlorocebus_sabaeus_excel_20200424/chlorochebus_sabaeus_w_cas9_scores_"

############### CAS OFF FINDER ################
CAS_OFF_EXCEL = "Chlorocebus_sabaeus_excel_20200424/Chlorocebus_sabaeus_filter_out_05_65_20200506/chlorochebus_sabaeus_w_cas9_scores_"
CAS_OFF_TXT = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_cas_NGG_off_"
CAS_OFF_TXT_TOP_20 = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_top20_cas_NGG_off_"
MAX_MISMATCH = 3
SUB_SET_NUM = 100
REF_SRV_PATH = "Input/FASTA/chlorocebus_sabaeus_chr"
PAM_N = "NNN"

# INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, SUB_SET_NUM, WORK_DIR + CAS_OFF_TXT, REF_SRV_PATH]
INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, 18, WORK_DIR + CAS_OFF_TXT_TOP_20, REF_SRV_PATH]

########### MERGE OFF_TAGET RESULT #############
CNT_RESULT_PATH = "CAS_OFF_FINDER/Output/CountResult/chlorochebus_sabaeus_top20_cas_NGG_off_"
# CNT_RESULT_PATH = "CAS_OFF_FINDER/Output/CountResult/ACE2_DPP4_TMPRSS2_"
CNT_RESULT_EXT = "_result_count.txt"
MAX_SEQ_NUM = 4
# GUIDE_SEQ_NUM = 20
FILE_CNT = 19
# FILE_CNT = 1
# FILE_CNT = SUB_SET_NUM

FILTER_OUT_OPT = [["1", "0", "0"], ["1", "0", ""], ["1", "", ""]]
CAS_SCORE_OPT = [50, 100]
FILTER_OUT_OPT_TOT = [FILTER_OUT_OPT, CAS_SCORE_OPT]

INITIAL_MRGE_OFF_TRGT = [WORK_DIR + CNT_RESULT_PATH, CNT_RESULT_EXT, FILE_CNT, PAM_N]
############### end setting env ################

def merge_off_target_total():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicsPrep()
    logic = Logic.Logics()

    off_trgt_dict = util.read_txt_to_dict(INITIAL_MRGE_OFF_TRGT)

    for i in range(1, 30):
        FILE_NAME_LIST.append(str(i))

    seq_cnt_group_by_crpt_id = {}
    re_off_trgt_dict = {}
    first_merge_dict = {}
    result_dict = {}
    for f_num in FILE_NAME_LIST:
        print("starting with file [" + str(f_num) + "]")
        df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, f_num)

        first_merge_dict = logic_prep.merge_excel_n_off_trgt(df_obj, off_trgt_dict)

        for cas_scr_opt in CAS_SCORE_OPT:
            for off_trg_opt_arr in FILTER_OUT_OPT:
                result_dict = logic.filter_out_by_rule(cas_scr_opt, off_trg_opt_arr, MAX_SEQ_NUM, first_merge_dict,
                                                       result_dict)

    seq_cnt_group_by_crpt_id = util.make_excel_w_off_trgt(WORK_DIR + CAS_OFF_EXCEL + "off_trgt_", "total", result_dict,
                                                          seq_cnt_group_by_crpt_id)

    re_off_trgt_dict = logic.get_re_off_target_seq(MAX_SEQ_NUM, seq_cnt_group_by_crpt_id, first_merge_dict,
                                                   re_off_trgt_dict)

    util.make_tab_txt_seq_cnt_group_by_crpt_id(WORK_DIR + CAS_OFF_EXCEL + "seq_cnt_group_by_crpt_id_",
                                               seq_cnt_group_by_crpt_id)

    util.make_tab_txt_re_off_target_seq(WORK_DIR + CAS_OFF_EXCEL + "re_off_target_seq_",
                                        re_off_trgt_dict)

def merge_off_target_indi():
    logic = Logic.Logics()

    for i in range(1, 30):
        FILE_NAME_LIST.append(str(i))

    logic.merge_off_target_indi(INITIAL_MRGE_OFF_TRGT, FILE_NAME_LIST, WORK_DIR, CAS_OFF_EXCEL, MAX_SEQ_NUM, FILTER_OUT_OPT_TOT)

def merge_deep_cas_9():
    util = Util.Utils()
    valid = Valid.Validations()

    for i in range(1, 30):
        FILE_NAME_LIST.append(str(i))

    for f_num in FILE_NAME_LIST:
        tmp_arr = util.read_cas9_scores_arr(WORK_DIR + RANK_FNAME, f_num)
        df_obj = util.read_excel_2_dataframe(WORK_DIR + EXCEL_FNAME, f_num)
        if valid.is_same_len(tmp_arr, df_obj, f_num):
            continue

        tmp_dict = util.merge_deep_cas_9_data(tmp_arr, df_obj)

        print(f_num + " is DONE")
        util.make_excel_srt_by_deepcas9_scr(WORK_DIR + SORTED_EXCEL_FNAME, f_num, tmp_dict)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
merge_off_target_total()
# merge_off_target_indi()
# merge_deep_cas_9()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))
