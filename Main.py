from time import clock

import Util
import Logic
############### start to set env ################
WORK_DIR = "D:/000_WORK/SongMyunJae_YuGooSang/20200417/WORK_DIR/"

SEQ_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome."

CDS_EACH_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chromosome."
CDS_FNAME = "Chlorocebus_sabaeus.ChlSab1.1.99.chr"

PAM_SEQ = ['NGG']
ADD_SEQ1_LEN = 4
SPACER_LEN = 20
ADD_SEQ2_LEN = 3
CLVG_AFTER_PAM = 3

FILE_NAME_LIST = ["X", "Y"]


INITIAL_SEQ = [PAM_SEQ, ADD_SEQ1_LEN, SPACER_LEN, ADD_SEQ2_LEN, CLVG_AFTER_PAM, WORK_DIR]
############### end setting env  ################

def main():
    util = Util.Utils()
    logic = Logic.Logics()

    for i in range(1,30):
        FILE_NAME_LIST.append(str(i))
    description_dict = util.read_dat_file_for_description(WORK_DIR + CDS_EACH_FNAME, FILE_NAME_LIST)

    cds_dict, sm_gene_diff_trscrpt = util.read_gtf_file_by_line_to_dict(WORK_DIR + CDS_FNAME, description_dict)

    logic.get_guide_ref(cds_dict, WORK_DIR + SEQ_FNAME, INITIAL_SEQ)

    # util.make_excel(result_dict, INITIAL_SEQ)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
main()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))
