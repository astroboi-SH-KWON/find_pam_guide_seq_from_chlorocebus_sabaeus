

class Validations:
    def __init__(self):
        self.tmp = ""

    def is_same_len(self, obj1, obj2, f_num):
        obj1_len = len(obj1)
        obj2_len = len(obj2)
        if obj1_len == obj2_len:
            print(f_num + " files are checked : [" + str(obj1_len) + "]")
            return False
        else:
            print(f_num + " is not matching, [" + str(obj1_len) + "] vs [" + str(obj2_len) + "]")
            return True
