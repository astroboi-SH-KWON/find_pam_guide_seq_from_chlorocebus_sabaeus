

import Util

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
        result_dict = {}
        idx = 1
        for key , vals in cds_dict.items():
            tmp_p_dict, tmp_m_dict = util.read_seq_dict(path, key, init)
            for val_dict in vals.values():
                if 'CDS' in val_dict:
                    cds_seq_arr = val_dict['CDS']
                    for cds_seq_exon_num in cds_seq_arr:
                        for i in range(cds_seq_exon_num[0] + 1, cds_seq_exon_num[1]):
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
                                idx = idx + 1
                            if i in tmp_m_dict:
                                result_dict[idx] = {}
                                result_dict[idx].update({'Target gene name': val_dict['Target gene name']})
                                result_dict[idx].update({'Ensembl transcript ID': val_dict['Ensembl transcript ID']})
                                result_dict[idx].update({'Ensembl Gene ID': val_dict['Ensembl Gene ID']})
                                if 'Description' in val_dict:
                                    result_dict[idx].update({'Description': val_dict['Description']})
                                result_dict[idx].update({'Position of Base After cut': i})
                                result_dict[idx].update({'Target context sequence': tmp_m_dict[i]})
                                result_dict[idx].update({'Strand': '-'})
                                result_dict[idx].update({'Exon Number': cds_seq_exon_num[2]})
                                idx = idx + 1

        return result_dict