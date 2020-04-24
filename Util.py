from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import re

import Logic

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"
        self.ext_xlsx = ".xlsx"

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

    """
    Chlorocebus_sabaeus.ChlSab1.1.99.chr.gtf
        this file doesn't have full gene name.
    :return
        {
        '8': {1: {'Target gene name': 'FBXO25', 'Description': 'F-box protein 25 [Source:HGNC', 'Ensembl transcript ID': 'ENSCSAT00000015163.1', 'Ensembl Gene ID': 'ENSCSAG00000017073.1', 'Strand': '+', 'CDS': ['189283 189416 3', '209061 209164 4', '210599 210648 5', '213558 213650 6', '228045 228138 7', '229300 229484 8', '236308 236490 9', '240929 241072 10', '245630 245656 11', '246617 246703 12']}
        , 2: {'Target gene name': 'ensembl', 'Ensembl transcript ID': 'ENSCSAT00000023666.1', 'Ensembl Gene ID': 'ENSCSAG00000023576.1', 'Strand': '+'}
        , 3: {'Target gene name': 'ERICH1', 'Description': 'glutamate rich 1 [Source:HGNC', 'Ensembl transcript ID': 'ENSCSAT00000015160.1', 'Ensembl Gene ID': 'ENSCSAG00000017066.1', 'Strand': '-', 'CDS': ['503913 503934 1', '488820 488966 2', '468987 469121 3', '452488 453222 4', '449438 449632 5', '445395 445465 6']}
        , 4: {'Target gene name': 'DLGAP2', 'Description': 'DLG associated protein 2 [Source:HGNC', 'Ensembl transcript ID': 'ENSCSAT00000015126.1', 'Ensembl Gene ID': 'ENSCSAG00000017056.1', 'Strand': '+', 'CDS': ['1255456 1256558 1', '1259920 1259929 2', '1263697 1263744 3', '1265494 1265508 4', '1269197 1269244 5', '1270408 1270443 6', '1284486 1284499 7', '1304565 1304712 8', '1309099 1309318 9', '1340355 1340704 10', '1341778 1341819 11', '1348440 1348525 12', '1350445 1350860 13', '1364929 1365020 14', '1370036 1370188 15', '1373314 1373529 16']}
        }
        , '7': {1: {'Target gene name': 'ensembl', 'Ensembl transcript ID': 'ENSCSAT00000003306.1', 'Ensembl Gene ID': 'ENSCSAG00000005272.1', 'Strand': '+', 'CDS': ['88663 89383 1', '90670 90818 2', '94479 94610 3', '94874 94961 4', '169533 169752 5', '255130 255406 6']}
        , 2: {'Target gene name': 'ensembl', 'Ensembl transcript ID': 'ENSCSAT00000003303.1', 'Ensembl Gene ID': 'ENSCSAG00000005268.1', 'Strand': '+', 'CDS': ['315966 316686 1', '318208 318356 2', '334463 334594 3', '338100 338187 4', '339022 339241 5', '343413 343689 6']}
        , 3: {'Target gene name': 'ensembl', 'Ensembl transcript ID': 'ENSCSAT00000003300.1', 'Ensembl Gene ID': 'ENSCSAG00000005266.1', 'Strand': '-'}
        , ...
        }
    """
    def read_gtf_file_by_line_to_dict(self, path, description_dict):
        tmp_dict = {}
        sm_gene_diff_trscrpt = []
        test_flag = 0
        with open(path + self.ext_gtf, "r") as f:
            idx = 0
            target_gene_name = ""
            while True:
                tmp_line = f.readline().replace("\n", "")

                if tmp_line != '':
                    if "#!" in tmp_line:
                        continue

                    line_arr = tmp_line.split("\t")

                    title = line_arr[2]
                    details_arr = line_arr[8].split(";")
                    chr_fname = line_arr[0]

                    if chr_fname == "MT":
                        continue

                    # # TODO delete when it's done
                    elif chr_fname in ["X","Y","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29"]:
                        continue
                    if chr_fname not in tmp_dict:
                        tmp_dict[chr_fname] = {}
                        idx = 0

                    if title == "gene":
                        target_gene_name = details_arr[2].split(" ")[2].replace('"', '')
                        test_flag = 0
                    elif title == "transcript":
                        idx = idx + 1
                        tmp_dict[chr_fname][idx] = {}
                        tmp_dict[chr_fname][idx].update({"Target gene name": target_gene_name})
                        if chr_fname in description_dict:
                            if target_gene_name in description_dict[chr_fname]:
                                tmp_dict[chr_fname][idx].update({"Description": description_dict[chr_fname][target_gene_name]})
                        ensembl_trscrpt_id = details_arr[2].split(" ")[2].replace('"', '') + "." + \
                                             details_arr[3].split(" ")[2].replace('"', '')
                        tmp_dict[chr_fname][idx].update({"Ensembl transcript ID": ensembl_trscrpt_id})
                        gene_id = details_arr[0].split(" ")[1].replace('"', '') + "." + details_arr[1].split(" ")[
                            2].replace('"', '')
                        tmp_dict[chr_fname][idx].update({"Ensembl Gene ID": gene_id})
                        tmp_dict[chr_fname][idx].update({"Strand": line_arr[6]})

                        test_flag = test_flag + 1
                        if test_flag > 1:
                            sm_gene_diff_trscrpt.append(
                                "chr [" + chr_fname + "] target_gene_nm : " + target_gene_name + ", trscrpt_id : " + ensembl_trscrpt_id + ", gene_id : " + gene_id)

                    elif title == "CDS":
                        if "CDS" in tmp_dict[chr_fname][idx]:
                            tmp_dict[chr_fname][idx]["CDS"].append(
                                line_arr[3] + " " + line_arr[4] + " " + details_arr[4].split(" ")[2].replace('"', ''))
                        else:
                            tmp_cds_seq_list = [
                                line_arr[3] + " " + line_arr[4] + " " + details_arr[4].split(" ")[2].replace('"', '')]
                            tmp_dict[chr_fname][idx].update({"CDS": tmp_cds_seq_list})
                else:
                    break

        return tmp_dict, sm_gene_diff_trscrpt

    """
    Chlorocebus_sabaeus.ChlSab1.1.99.chromosome.1.dat, 2.dat, 3.dat ... Y.dat

    :return
        {'X': {'DHRSX': 'dehydrogenase/reductase X-linked [Source:HGNC'
                , 'ZBED1': 'zinc finger BED-type containing 1 [Source:NCBI'
                , 'XG': 'Xg glycoprotein (Xg blood group) [Source:NCBI', 'SRY': 'sex determining region Y [Source:NCBI'
                }
        , '2': {'PCMTD2': 'protein-L-isoaspartate (D-aspartate) O-'
                , 'MYT1': 'myelin transcription factor 1 [Source:HGNC'
                ,  ...}
    """
    def read_dat_file_for_description(self, path, file_num_arr):
        tmp_dict = {}
        for i in file_num_arr:
            tmp_dict[i] = {}
            with open(path + str(i) + self.ext_dat, "r") as f:

                while True:
                    tmp_line = f.readline().replace("\n", "")

                    if tmp_line != '':
                        # /locus_tag= is Target gene name, /note= right under /locus_tag= is Description of Target gene name
                        if "/locus_tag=" in tmp_line:
                            target_gene_name = tmp_line.replace(" ", "").replace("/locus_tag=", "").replace('"', '')
                            tmp_description = f.readline().replace("\n", "")
                            if "/note=" in tmp_description:
                                tmp_description = tmp_description.replace("/note=","").replace('"', '').lstrip()
                                tmp_dict[i].update({target_gene_name: tmp_description})
                            else:
                                print(target_gene_name + " in file [" + str(i) + "] doesn't have full description ")
                    else:
                        break
        return tmp_dict


    """
    Chlorocebus_sabaeus.ChlSab1.1.99.chromosome.1.dat, 2.dat, 3.dat ... Y.dat
        these files don't have exon number, and their CDS sequences are odd
    """
    def read_dat_file_by_line_to_dict(self, path):
        tmp_dict = {}
        with open(path + self.ext_dat, "r") as f:
            idx = 0
            flag = True
            while True:
                tmp_line = f.readline().replace("\n", "")

                while flag:
                    comment_line = f.readline().replace("\n", "")
                    if "FEATURES" in comment_line:
                        flag = False
                        break

                if tmp_line != '':
                    if "ORIGIN" in tmp_line:
                        break
                    # if "gene=ENSCSAG00000023543.1" in tmp_line:
                    #     break
                    if "    gene " in tmp_line:
                        idx = idx + 1
                        tmp_dict[idx] = {}
                        # print("tmp_line : " + tmp_line)
                        target_str_arr = tmp_line.split("..")
                        tmp_num0 = re.sub("\D", "", target_str_arr[0])
                        tmp_num1 = re.sub("\D", "", target_str_arr[1])
                        if tmp_num0 > tmp_num1:
                            tmp_dict[idx].update({"Strand": "+"})
                            print(str(idx) + " it's plus")
                        else:
                            tmp_dict[idx].update({"Strand": "-"})

                    if "/gene=" in tmp_line:
                        if '"' not in tmp_line:
                            tmp_dict[idx].update({"Ensembl Gene ID": tmp_line.replace(" ", "").replace("/gene=", "")})
                    if "/locus_tag=" in tmp_line:
                        tmp_dict[idx].update({"Target gene name": tmp_line.replace(" ", "").replace("/locus_tag=", "").replace('"', '')})
                    if "/note=" in tmp_line:
                        if "Description" not in tmp_dict[idx]:
                            tmp_dict[idx].update({"Description": tmp_line.replace("/note=","").replace('"', '').lstrip()})
                    if "/standard_name=" in tmp_line:
                        tmp_dict[idx].update({"Ensembl transcript ID": tmp_line.replace(" ", "").replace("/standard_name=", "").replace('"', '')})
                        cds_line = ""
                        while True:
                            cds_tmp_line = f.readline().replace("\n", "")
                            # case of no CDS array
                            if "    gene " in cds_tmp_line:
                                idx = idx + 1
                                tmp_dict[idx] = {}
                                # print("cds_tmp_line : " + cds_tmp_line)
                                target_str_arr = cds_tmp_line.split("..")
                                tmp_num0 = re.sub("\D", "", target_str_arr[0])
                                tmp_num1 = re.sub("\D", "", target_str_arr[1])
                                if tmp_num0 > tmp_num1:
                                    tmp_dict[idx].update({"Strand": "+"})
                                    print(str(idx) + " it's plus")
                                else:
                                    tmp_dict[idx].update({"Strand": "-"})
                                break

                            cds_tmp_line = cds_tmp_line.replace(" ", "").replace("CDSjoin(", "").replace("CDS", "").replace("complement(", "").replace(")", "")
                            if "gene" in cds_tmp_line:
                                cds_line_arr = cds_line.split(",")
                                tmp_dict[idx].update({"CDS": cds_line_arr})
                                break
                            if "protein_id" in cds_tmp_line:
                                cds_line_arr = cds_line.split(",")
                                tmp_dict[idx].update({"CDS": cds_line_arr})
                                break
                            if "note" in cds_tmp_line:
                                cds_line_arr = cds_line.split(",")
                                tmp_dict[idx].update({"CDS": cds_line_arr})
                                break
                            if "translation" in cds_tmp_line:
                                cds_line_arr = cds_line.split(",")
                                tmp_dict[idx].update({"CDS": cds_line_arr})
                                break
                            cds_line = cds_line + cds_tmp_line

                else:
                    break
        return tmp_dict

    def read_seq_dict(self, path, fil_num_str, init_arr):
        pam_arr = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        clvg_after_pam = init_arr[4]
        tmp_p_dict = {}
        tmp_m_dict = {}
        logic = Logic.Logics()

        for pam_str in pam_arr:
            pam_len = len(pam_str)
            std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len
            with open(path + fil_num_str + self.ext_fa, "r") as f:
                header = f.readline()  # header ignored : >chr19
                # print("header : " + header)

                idx = 1
                tmp_p_str = ""
                tmp_m_str = ""

                while True:
                    c = f.read(1)
                    if c == "":
                        break
                    if "\n" in c:
                        continue
                    elif "\r" in c:
                        continue
                    # # TODO delete when it's done
                    # if idx > 6425128:
                    #     break

                    tmp_p_str = tmp_p_str + c.upper()
                    tmp_m_str = tmp_m_str + logic.get_complementary(c.upper())

                    if len(tmp_p_str) > std_tot_len:
                        tmp_p_str = tmp_p_str[-std_tot_len:]
                        tmp_m_str = tmp_m_str[-std_tot_len:]

                    if len(tmp_p_str) == std_tot_len:
                        # filter out 'N' seq
                        if 'N' not in tmp_p_str:
                            if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_str):
                                clvg_p_idx = idx - (pam_len + add_seq2_len + clvg_after_pam - 1)
                                tmp_p_dict[clvg_p_idx] = tmp_p_str
                                # print("tmp_p_str : " + tmp_p_str + ", clvg_p_idx : " + str(clvg_p_idx))

                            if logic.match(0, tmp_m_str[add_seq2_len:add_seq2_len + pam_len], pam_str[::-1]):
                                clvg_m_idx = idx - (std_tot_len - 1) + add_seq2_len + pam_len + clvg_after_pam
                                # tmp_m_dict[clvg_m_idx] = tmp_m_str
                                # expression is + strand type even in - strand
                                tmp_m_dict[clvg_m_idx] = tmp_p_str + " " +tmp_m_str
                                # print("tmp_m_str : " + tmp_m_str + ", clvg_m_idx : " + str(clvg_m_idx))

                    idx = idx + 1

        return tmp_p_dict, tmp_m_dict

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
            sheet.cell(row=row, column=1, value=key)
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


