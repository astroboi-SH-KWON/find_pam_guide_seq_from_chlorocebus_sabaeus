from time import clock

import Util
import Logic
############### start to set env ################
WORK_DIR = "D:/000_WORK/SongMyunJae_YuGooSang/20200417/WORK_DIR/"

FA_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome."
FA_FNAME_EXT = ".fa"

CDS_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chromosome."
CDS_FNAME_EXT = ".dat"


PAM_SEQ = ['NGG']



INITIAL_ANALYSIS = [PAM_SEQ, WORK_DIR + RAW_DATA_DIR, MAX_MIS_CNT]
MAKE_EXCEL = [WORK_DIR, COL_NAME, MAX_MIS_CNT]
############### end setting env  ################

