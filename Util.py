from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import random
import math

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def read_cas9_scores_arr(self, path, f_num):
        tmp_str = ""
        with open(path + f_num + self.ext_txt, "r") as f:
            line_001 = f.readline().replace("\n", "")
            line_002 = f.readline().replace("\n", "")
            line_003 = f.readline().replace("\n", "")
            line_004 = f.readline().replace("\n", "")
            tmp_str = f.readline().replace("\n", "")
        print(f_num + " 's DeepCas9 score is read")
        return tmp_str.replace("(", "").replace(")", "").split(",")

    def read_excel_2_dataframe(self, path, f_name):
        return pd.read_excel(path + f_name + self.ext_xlsx)

    def read_dat_file_by_line_to_list(self,path):
        tmp_list = []
        with open(path + self.ext_dat, "r") as f:
            gene_str = ""
            CDC_str = ""
            idx = 0
            while True:
                tmp_line = f.readline().replace("\n", "")

                if tmp_line != '':
                    if "gene=ENSCSAG00000016985.1" in tmp_line:
                        break

                    if " gene " in tmp_line:
                        idx = idx + 1
                        print(str(idx) + ":::::::::::::::::::::::::::::::::::")

                    # print(tmp_line)
                    tmp_list.append(tmp_line)
                else:
                    break
        return tmp_list

    def make_txt_with_list(self, path, data_list):
        tmp_set = set()
        with open(path + "_" + str(clock()) + self.ext_txt, 'a') as f:
            for data_str in data_list:
                f.write(data_str + "\n")
                gene_id = data_str.split(" ")[-1]
                tmp_set.add(gene_id)
        return tmp_set

    def make_excel(self, result_dict, init, f_num):
        pam_len = len(init[0][0])
        add_1_len = init[1]
        spcr_len = init[2]
        add_2_len = init[3]
        clvg_after_pam = init[4]
        path = init[5]

        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        sheet.cell(row=row, column=1, value="INDEX")
        sheet.cell(row=row, column=2, value='Target gene name')
        sheet.cell(row=row, column=3, value='Description')
        sheet.cell(row=row, column=4, value='Ensembl transcript ID')
        sheet.cell(row=row, column=5, value='Ensembl Gene ID')
        sheet.cell(row=row, column=6, value='Position of Base After cut')
        sheet.cell(row=row, column=7, value='Strand')
        sheet.cell(row=row, column=8, value='sgRNA Target sequence')
        sheet.cell(row=row, column=9, value='Target context sequence')
        sheet.cell(row=row, column=10, value='PAM')
        sheet.cell(row=row, column=11, value='order sgRNA Target sequence')
        sheet.cell(row=row, column=12, value='order Target context sequence')
        sheet.cell(row=row, column=13, value='order PAM')
        sheet.cell(row=row, column=14, value='Exon Number')
        sheet.cell(row=row, column=15, value='DeepCas9 score')
        sheet.cell(row=row, column=16, value='target site (cleavage site)')

        for key, val_arr in result_dict.items():
            row = row + 1
            sheet.cell(row=row, column=1, value=(row - 1))
            if "ensembl".upper() != val_arr['Target gene name'].upper():
                sheet.cell(row=row, column=2, value=val_arr['Target gene name'])
            if 'Description' in val_arr:
                sheet.cell(row=row, column=3, value=val_arr['Description'])
            sheet.cell(row=row, column=4, value=val_arr['Ensembl transcript ID'])
            sheet.cell(row=row, column=5, value=val_arr['Ensembl Gene ID'])
            sheet.cell(row=row, column=6, value=val_arr['Position of Base After cut'])
            strand = val_arr['Strand']
            sheet.cell(row=row, column=7, value=strand)
            seq_str = val_arr['Target context sequence']
            if strand == '+':
                sheet.cell(row=row, column=8, value=seq_str[add_1_len:add_1_len + spcr_len])
                sheet.cell(row=row, column=9, value=seq_str)
                sheet.cell(row=row, column=10, value=seq_str[add_1_len + spcr_len:-add_2_len])
                sheet.cell(row=row, column=11, value=seq_str[add_1_len:add_1_len + spcr_len])
                sheet.cell(row=row, column=12, value=seq_str)
                sheet.cell(row=row, column=13, value=seq_str[add_1_len + spcr_len:-add_2_len])
            else:
                comp_seq_str = val_arr['Target context anti sequence']
                sheet.cell(row=row, column=8, value=seq_str[add_2_len + pam_len:add_2_len + pam_len + spcr_len])
                sheet.cell(row=row, column=9, value=seq_str)
                sheet.cell(row=row, column=10, value=seq_str[add_2_len:add_2_len + pam_len])
                sheet.cell(row=row, column=11, value=comp_seq_str[add_2_len + pam_len:add_2_len + pam_len + spcr_len][::-1])
                sheet.cell(row=row, column=12, value=comp_seq_str[::-1])
                sheet.cell(row=row, column=13, value=comp_seq_str[add_2_len:add_2_len + pam_len][::-1])
            sheet.cell(row=row, column=14, value=val_arr['Exon Number'])
            sheet.cell(row=row, column=16, value=val_arr['Ratio'])
        workbook.save(filename=path + "chlorochebus_sabaeus_" + f_num + "_" + str(clock()) + self.ext_xlsx)

    def make_excel_srt_by_deepcas9_scr(self, path, f_num, tmp_dict):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        for key, vals in tmp_dict.items():
            for tmp_list in sorted(vals, key=lambda tmp_list: tmp_list[0],reverse=True):
                col = 14
                for tmp_val in tmp_list:
                    if col >= 14:
                        sheet.cell(row=row, column=col, value=tmp_val)
                        if col >= 15:
                            col = 0
                    else:
                        if col < 3:
                            sheet.cell(row=row, column=col, value=tmp_val)
                        elif col == 3:
                            sheet.cell(row=row, column=col, value=key)
                            sheet.cell(row=row, column=col + 1, value=tmp_val)
                        else:
                            sheet.cell(row=row, column=col + 1, value=tmp_val)
                    col = col + 1
                row = row + 1

        # workbook.save(filename=path + f_num + "_" + str(clock()) + self.ext_xlsx)
        workbook.save(filename=path + f_num + self.ext_xlsx)

    def merge_deep_cas_9_data(self, tmp_arr, df_obj):
        tmp_dict = {'Ensembl transcript ID': [
            ['DeepCas9 score', 'target site (cleavage site)', 'Target gene name', 'Description', 'Ensembl Gene ID',
             'Position of Base After cut', 'Strand', 'sgRNA Target sequence', 'Target context sequence', 'PAM',
             'order sgRNA Target sequence', 'order Target context sequence', 'order PAM', 'Exon Number']]}
        data_obj = df_obj.to_dict()
        for i in range(len(tmp_arr)):
            # filtering out if target site (cleavage site) < 0.05, > 0.65
            if data_obj['target site (cleavage site)'][i] < 0.05:
                continue
            if data_obj['target site (cleavage site)'][i] > 0.65:
                continue
            trnscrpt_id = data_obj['Ensembl transcript ID'][i]
            tmp_list = []
            tmp_list.append(float(tmp_arr[i]))
            tmp_list.append(data_obj['target site (cleavage site)'][i])
            tmp_list.append(data_obj['Target gene name'][i])
            tmp_list.append(data_obj['Description'][i])
            tmp_list.append(data_obj['Ensembl Gene ID'][i])
            tmp_list.append(data_obj['Position of Base After cut'][i])
            tmp_list.append(data_obj['Strand'][i])
            tmp_list.append(data_obj['sgRNA Target sequence'][i])
            tmp_list.append(data_obj['Target context sequence'][i])
            tmp_list.append(data_obj['PAM'][i])
            tmp_list.append(data_obj['order sgRNA Target sequence'][i])
            tmp_list.append(data_obj['order Target context sequence'][i])
            tmp_list.append(data_obj['order PAM'][i])
            tmp_list.append(data_obj['Exon Number'][i])
            if trnscrpt_id not in tmp_dict.keys():
                tmp_dict[trnscrpt_id] = [tmp_list]
            else:
                tmp_dict[trnscrpt_id].append(tmp_list)
        return tmp_dict

    """
    make_txt_cas_off_finder_input : .txt input file for CAS_OFF_FINDER
    :param
        path
        seq_set_pool
        mis_mtch_num : maximum # of mismatch count
        sub_set_num : sub set # of seq_set_pool to make .txt input files separately
    """
    def make_txt_cas_off_finder_input(self, seq_set_pool, init):
        pam_seq = init[0]
        spacr_len = init[1]
        mis_mtch_num = init[2]
        sub_set_num = init[3]
        path = init[4]
        ref_srv_path = init[5]

        raw_each_sub_len = len(seq_set_pool) / sub_set_num
        each_sub_len = int(raw_each_sub_len)
        # * 1000) / 1000 : raw_each_sub_len 값이 반욜림되는 경우의 error 방지
        raw_last_sub_len_to_add = int(((raw_each_sub_len - each_sub_len) * sub_set_num) * 1000)/1000
        last_sub_len_to_add = math.ceil(raw_last_sub_len_to_add)
        print("total len : " + str(len(seq_set_pool)) + ", raw_each_sub_len : " + str(raw_each_sub_len))
        print("raw_last_sub_len_to_add : " + str(raw_last_sub_len_to_add) + " , last_sub_len_to_add : " + str(last_sub_len_to_add))

        for fname_idx in range(sub_set_num):
            if fname_idx == sub_set_num -1:
                each_sub_len += last_sub_len_to_add
            tmp_sub_set = set(random.sample(seq_set_pool, each_sub_len))
            seq_set_pool -= tmp_sub_set
            print("tmp_sub_set_" + str(fname_idx) + " len : " + str(len(tmp_sub_set)))
            print("rest total len : " + str(len(seq_set_pool)))
            print(" ")
            with open(path + str(fname_idx) + self.ext_txt, 'a') as f:
                f.write(ref_srv_path + "\n")
                f.write("N"*spacr_len + pam_seq + "\n")
                for data_str in tmp_sub_set:
                    f.write(data_str + " " + str(mis_mtch_num) +"\n")

    def read_txt_to_dict(self, path, file_ext, cnt):
        result_dict = {}
        for idx in range(cnt):
            with open(path + str(idx) + file_ext, "r") as f:
                header = f.readline()
                # print("header : " + header)# Sequence	0	1	2	3

                while True:
                    read_line = f.readline().replace("\n","")
                    if read_line == "":
                        break
                    read_line_arr = read_line.split("\t")
                    tmp_list = []
                    for i in range(1, len(read_line_arr)):
                        tmp_list.append(int(read_line_arr[i]))
                    result_dict[read_line_arr[0]] = tmp_list

        return result_dict

    def make_excel_off_target_data(self, off_trgt_dict, df_obj, path, f_num):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        sheet.cell(row=row, column=1, value="INDEX")
        sheet.cell(row=row, column=2, value='Target gene name')
        sheet.cell(row=row, column=3, value='Description')
        sheet.cell(row=row, column=4, value='Ensembl transcript ID')
        sheet.cell(row=row, column=5, value='Ensembl Gene ID')
        sheet.cell(row=row, column=6, value='Position of Base After cut')
        sheet.cell(row=row, column=7, value='Strand')
        sheet.cell(row=row, column=8, value='sgRNA Target sequence')
        sheet.cell(row=row, column=9, value='Target context sequence')
        sheet.cell(row=row, column=10, value='PAM')
        sheet.cell(row=row, column=11, value='order sgRNA Target sequence')
        sheet.cell(row=row, column=12, value='order Target context sequence')
        sheet.cell(row=row, column=13, value='order PAM')
        sheet.cell(row=row, column=14, value='Exon Number')
        sheet.cell(row=row, column=15, value='DeepCas9 score')
        sheet.cell(row=row, column=16, value='target site (cleavage site)')
        col_head = 17
        for mis_match_cnt in range(len(off_trgt_dict[list(off_trgt_dict)[0]])):
            sheet.cell(row=row, column=col_head, value=str(mis_match_cnt))
            col_head += 1

        row += 1
        data_obj = df_obj.to_dict()
        idx = 0
        for i in range(len(data_obj['Ensembl transcript ID'])):
            ordr_trgt_seq = data_obj['order sgRNA Target sequence'][i]
            ordr_pam = data_obj['order PAM'][i]

            if ordr_trgt_seq + ordr_pam in off_trgt_dict:
                off_trgt_arr = off_trgt_dict[ordr_trgt_seq + ordr_pam]
                # filtering out if result of (mis_match = 0) > 1
                if off_trgt_arr[0] > 1:
                    continue
                col = 17
                for result_cnt in off_trgt_arr:
                    sheet.cell(row=row, column=col, value=str(result_cnt))
                    col += 1
            idx += 1
            sheet.cell(row=row, column=1, value=str(idx))
            sheet.cell(row=row, column=2, value=data_obj['Target gene name'][i])
            sheet.cell(row=row, column=3, value=data_obj['Description'][i])
            sheet.cell(row=row, column=4, value=data_obj['Ensembl transcript ID'][i])
            sheet.cell(row=row, column=5, value=data_obj['Ensembl Gene ID'][i])
            sheet.cell(row=row, column=6, value=data_obj['Position of Base After cut'][i])
            sheet.cell(row=row, column=7, value=data_obj['Strand'][i])
            sheet.cell(row=row, column=8, value=data_obj['sgRNA Target sequence'][i])
            sheet.cell(row=row, column=9, value=data_obj['Target context sequence'][i])
            sheet.cell(row=row, column=10, value=data_obj['PAM'][i])
            sheet.cell(row=row, column=11, value=ordr_trgt_seq)
            sheet.cell(row=row, column=12, value=data_obj['order Target context sequence'][i])
            sheet.cell(row=row, column=13, value=ordr_pam)
            sheet.cell(row=row, column=14, value=data_obj['Exon Number'][i])
            sheet.cell(row=row, column=15, value=data_obj['DeepCas9 score'][i])
            sheet.cell(row=row, column=16, value=data_obj['target site (cleavage site)'][i])

            row += 1

        workbook.save(filename=path + f_num + "_" + str(clock()) + self.ext_xlsx)
