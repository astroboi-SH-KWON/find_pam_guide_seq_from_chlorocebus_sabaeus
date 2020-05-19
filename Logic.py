

import Util
import LogicPrep

class Logics:
    def __init__(self):
        self.tmp = ""

    def get_complementary(self, c):
        if c == 'C':
            return "G"
        elif c == 'A':
            return "T"
        elif c == 'T':
            return "A"
        elif c == 'G':
            return "C"
        elif c == 'N':
            return "N"
        else:
            print("get_complementary ERROR .... char is [" + c + "]")
            exit()

    """
    match : match sequence with same length strings
    :param
        i : index of seq
        dna_seq : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """
    def match(self,i, dna_seq, rule_str):
        if len(dna_seq) == i:
            return True
        if self.checkSeqByChar(dna_seq[i], rule_str[i]):
            return self.match(i + 1, dna_seq, rule_str)
        else:
            return False

    """
    checkSeqByChar : match sequences by char with rules
    :param
        dna_char :
        rule_char : rules with "A", "C", "G", "T", "U", "N", "R",...
    :return
        boolean
    """
    def checkSeqByChar(self,dna_char, rule_char):
        flag = False
        if rule_char == 'N':
            return True
        elif rule_char in 'ACGTU':
            if dna_char == rule_char:
                return True
        elif rule_char == 'R':
            if dna_char in 'AG':
                return True
        # elif rule_char == 'r':
        #     if dna_char in 'CT':
        #         return True
        """
        add more rules of "ACTGU"
        """
        return flag

    def get_guide_ref(self, cds_dict, path, init):
        util = Util.Utils()
        logic_pre = LogicPrep.LogicsPrep()
        # result_dict = {}
        idx = 1
        for key , vals in cds_dict.items():

            result_dict = {}
            tmp_p_dict, tmp_m_dict = logic_pre.read_seq_dict(path, key, init)
            for val_dict in vals.values():
                if 'CDS' in val_dict:
                    cds_seq_arr = val_dict['CDS']

                    total_cds_len = 0
                    for cds_tmp in cds_seq_arr:
                        cds_seq_tmp = cds_tmp.split(" ")
                        total_cds_len = total_cds_len + (int(cds_seq_tmp[1]) - int(cds_seq_tmp[0]) + 1)

                    prent_cds_len = 0
                    for cds_seq in cds_seq_arr:
                        cds_seq_exon_num = cds_seq.split(" ")
                        prent_cds_len = prent_cds_len + (int(cds_seq_exon_num[1]) - int(cds_seq_exon_num[0]) + 1)
                        for i in range(int(cds_seq_exon_num[0]) + 1, int(cds_seq_exon_num[1]) + 1):
                            if i in tmp_p_dict:
                                result_dict[idx] = {}
                                result_dict[idx].update({'Target gene name': val_dict['Target gene name']})
                                result_dict[idx].update({'Ensembl transcript ID': val_dict['Ensembl transcript ID']})
                                result_dict[idx].update({'Ensembl Gene ID': val_dict['Ensembl Gene ID']})
                                if 'Description' in val_dict:
                                    result_dict[idx].update({'Description': val_dict['Description']})
                                result_dict[idx].update({'Position of Base After cut': i})
                                result_dict[idx].update({'Target context sequence': tmp_p_dict[i]})
                                result_dict[idx].update({'Strand': '+'})
                                result_dict[idx].update({'Exon Number': cds_seq_exon_num[2]})
                                # the pos ratio out of total CDS depends on which strand gene occurs
                                if '+' == val_dict['Strand']:
                                    this_len = prent_cds_len - (int(cds_seq_exon_num[1]) - i + 1)
                                    result_dict[idx].update({'Ratio': this_len / total_cds_len})
                                else:
                                    this_len = prent_cds_len - (i - int(cds_seq_exon_num[0]) + 1)
                                    result_dict[idx].update({'Ratio': this_len / total_cds_len})
                                idx = idx + 1
                            if i in tmp_m_dict:
                                result_dict[idx] = {}
                                result_dict[idx].update({'Target gene name': val_dict['Target gene name']})
                                result_dict[idx].update({'Ensembl transcript ID': val_dict['Ensembl transcript ID']})
                                result_dict[idx].update({'Ensembl Gene ID': val_dict['Ensembl Gene ID']})
                                if 'Description' in val_dict:
                                    result_dict[idx].update({'Description': val_dict['Description']})
                                result_dict[idx].update({'Position of Base After cut': i})
                                result_dict[idx].update({'Target context sequence': tmp_m_dict[i].split(" ")[0]})
                                result_dict[idx].update({'Target context anti sequence': tmp_m_dict[i].split(" ")[1]})
                                result_dict[idx].update({'Strand': '-'})
                                result_dict[idx].update({'Exon Number': cds_seq_exon_num[2]})
                                # the pos ratio out of total CDS depends on which strand gene occurs
                                if '+' == val_dict['Strand']:
                                    this_len = prent_cds_len - (int(cds_seq_exon_num[1]) - i + 1)
                                    result_dict[idx].update({'Ratio': this_len / total_cds_len})
                                else:
                                    this_len = prent_cds_len - (i - int(cds_seq_exon_num[0]) + 1)
                                    result_dict[idx].update({'Ratio': this_len / total_cds_len})
                                idx = idx + 1
            print("DONE file_" + key)
            util.make_excel(result_dict, init, key)
        # return result_dict

    """
    Guide 선정 기준
    * Mismatch 3개인 것의 개수는 고려하지 않습니다. (너무 많기 때문)
    1순위 - mismatch 1/2 개수 0인 것에서 Cas9 score 50점 이상
    2순위 -  Mismatch 1 개수 0인 것에서 Cas9 score 50점 이상 (mismatch 2bp까지 허용)
    3순위 - Mismatch 고려 안하고 Cas9 score 50점 이상 (mismatch 1bp까지 허용)

    여기서 4개 채워지지 않은 gene (정확히는 transcript)들은 Cas9 score 
        상위 200개로 다시 list up해서 off-target 돌려서 4개 guide 채우기
    
    :param
        cas9_scre_litmit : min of DeepCas9 score. opt for filtering out
        off_trg_opt_arr : off_target option array ex) ["1", "0", ""]
    """
    def filter_out_by_rule(self, cas9_scre_litmit, off_trg_opt_arr, max_num, data_dict, result_dict):

        for trnscrpt_id, list_arr in data_dict.items():

            is_not_first_id = True
            if trnscrpt_id in result_dict:
                if len(result_dict[trnscrpt_id]) > max_num - 1:
                    continue
            else:
                result_dict.update({trnscrpt_id: []})
                is_not_first_id = False

            for val_arr in list_arr:
                if len(result_dict[trnscrpt_id]) > max_num - 1:
                    break

                deep_cas9_score = val_arr[12]
                off_trgt_arr = val_arr[14]

                if len(off_trgt_arr) == 0:
                    continue

                if cas9_scre_litmit == 100:
                    pass
                elif deep_cas9_score < cas9_scre_litmit:
                    continue

                off_trgt_flag = False
                for idx in range(len(off_trg_opt_arr)):
                    off_cnt_str = off_trg_opt_arr[idx]
                    if off_cnt_str != "":
                        if off_trgt_arr[idx] > int(off_cnt_str):
                            off_trgt_flag = True
                            break

                if off_trgt_flag:
                    continue
                if is_not_first_id:
                    exist_flag = False
                    for tmp_arr in result_dict[trnscrpt_id]:
                        if tmp_arr == val_arr:
                            exist_flag = True
                            break
                    if not exist_flag:
                        result_dict[trnscrpt_id].append(val_arr)
                else:
                    result_dict[trnscrpt_id].append(val_arr)

        return result_dict

    """
    :param
        seq_cnt_group_by_crpt_id = {trnscrpt_id: # of seq}
    """
    def get_re_off_target_seq(self, max_num, seq_cnt_group_by_crpt_id, data_dict, result_dict):
        for trnscrpt_id, arr_list in data_dict.items():
            if seq_cnt_group_by_crpt_id[trnscrpt_id] < max_num:
                for seq_arr in arr_list:
                    if trnscrpt_id in result_dict:
                        result_dict[trnscrpt_id].append(seq_arr[8])
                    else:
                        result_dict.update({trnscrpt_id: [seq_arr[8]]})
        return result_dict

    def merge_off_target_indi(self, init_off_trgt, f_nm_list, wrk_dir, excel_path, max_num, filt_out_opt):
        util = Util.Utils()
        logic_prep = LogicPrep.LogicsPrep()

        off_trgt_opt = filt_out_opt[0]
        cas_scr_opt = filt_out_opt[1]

        off_trgt_dict = util.read_txt_to_dict(init_off_trgt)

        seq_cnt_group_by_crpt_id = {}
        re_off_trgt_dict = {}
        for f_num in f_nm_list:
            print("starting with file [" + str(f_num) + "]")
            df_obj = util.read_excel_2_dataframe(wrk_dir + excel_path, f_num)

            first_merge_dict = logic_prep.merge_excel_n_off_trgt(df_obj, off_trgt_dict)

            result_dict = {}
            for cas_opt in cas_scr_opt:
                for off_trg_opt_arr in off_trgt_opt:
                    result_dict = self.filter_out_by_rule(cas_opt, off_trg_opt_arr, max_num, first_merge_dict, result_dict)

            seq_cnt_group_by_crpt_id = util.make_excel_w_off_trgt(wrk_dir + excel_path + "off_trgt_", f_num,
                                                                  result_dict, seq_cnt_group_by_crpt_id)

            re_off_trgt_dict = self.get_re_off_target_seq(max_num, seq_cnt_group_by_crpt_id, first_merge_dict,
                                                          re_off_trgt_dict)
            print("done with file [" + str(f_num) + "]\n")

        util.make_tab_txt_seq_cnt_group_by_crpt_id(wrk_dir + excel_path + "seq_cnt_group_by_crpt_id_",
                                                   seq_cnt_group_by_crpt_id)

        util.make_tab_txt_re_off_target_seq(wrk_dir + excel_path + "re_off_target_seq_",
                                            re_off_trgt_dict)