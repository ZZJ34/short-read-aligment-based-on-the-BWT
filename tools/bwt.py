
class Suffix:
    def __init__(self, text, position):
        self.text = text
        self.pos = position

def text_key(a): 
    return a.text

def reference_array2str(reference_array):
        reference_str = ''

        # 按照数组的顺序结尾添加结束符组成字符串
        for reference_item in reference_array:
            reference_str = reference_str + "%s$" % reference_item.upper()

        return reference_str




# BWA 变换生成的 FM-index 数据
class BWA_FM_index:
    # BWT B, the occurrence array O(a,i),  suffix array S(i)
    def __init__(self, reference_array):

        reference = reference_array2str(reference_array)
        
        # 创建各种数据的存储类型
        rotation_list, suffix_array, bwt, = [list() for i in range(3)]
        C, Occ= [dict() for i in range(2)]

        # 创建一个无序不重复元素集
        alphabet = { '$', 'A', 'C', 'G', 'T', 'X'} 

        # 初始化 C, Occ (C, Occ 是 dict) 
        for char in alphabet:
            C[char] = 0  
            Occ[char] = list()  

        # 生成基因组参考序列的循环移位序列，并存储于 rotation_list
        for i in range(len(reference)):
            new_rotation = "%s%s" % (reference[i:], reference[0:i])
            struct = Suffix(new_rotation, i)

            rotation_list.append(struct)

        # 生成C 统计结束符
        for symbol in alphabet:
            for i in range(len(reference)):
                if reference[i] < symbol:
                    C[symbol] = C[symbol] + 1

        # 将循环移位序列根据 text 的字典序排序 
        # 注意这里先不考虑 N, 可以将 N 替换成 X 再进行排序
        rotation_list.sort(key=text_key)
        
        # for item in rotation_list:
        #     print(item.text)

        # 生成 suffix array 和 BWT
        for item in rotation_list:
            suffix_array.append(item.pos)
            bwt.append(item.text[-1:])

        # 生成 reference 的 Occ
        for item in rotation_list:
            for symbol in alphabet:
                if len(Occ[symbol]) == 0:
                    prev = 0
                else:
                    prev = Occ[symbol][-1]

                if item.text[-1:] == symbol:
                    Occ[symbol].append(prev + 1)
                else:
                    Occ[symbol].append(prev)

        self.data = dict()
        self.data['S'] = suffix_array
        self.data['B'] = bwt
        self.data['C'] = C
        self.data['O'] = Occ
        self.text = reference
        self.text_num = len(reference_array) 

    #注意：Occ(a,-1) = 0
    def Occ(self, char, index):
        if index < 0:
            return 0
        else:
            return self.data['O'][char][index]



# BWA 变换生成的 FMD-index 数据
class BWA_FMD_index:
    # BWT B, the occurrence array O(a,i),  suffix array S(i)
    def __init__(self, reference_array):
        
        # 互补匹配 dict
        self.reverse_dict={
            "A" : "T",
            "T" : "A",
            "C" : "G",
            "G" : "C",
            "X" : "X"
        }

        reference = reference_array2str(self.add_reverse_reference(reference_array))

        # 创建各种数据的存储类型
        rotation_list, suffix_array, bwt, = [list() for i in range(3)]
        C, Occ= [dict() for i in range(2)]

        # 创建一个无序不重复元素集
        alphabet = { '$', 'A', 'C', 'G', 'T', 'X'} 

        # 初始化 C, Occ (C, Occ 是 dict) 
        for char in alphabet:
            C[char] = 0  
            Occ[char] = list()  

        # 生成基因组参考序列的循环移位序列，并存储于 rotation_list
        for i in range(len(reference)):
            new_rotation = "%s%s" % (reference[i:], reference[0:i])
            struct = Suffix(new_rotation, i)
            rotation_list.append(struct)

        # 生成C 统计结束符
        for symbol in alphabet:
            for i in range(len(reference)):
                if reference[i] < symbol:
                    C[symbol] = C[symbol] + 1

        # 将循环移位序列根据 text 的字典序排序 
        # 注意这里先不考虑 N, 可以将 N 替换成 X 再进行排序
        rotation_list.sort(key=text_key)
        
        # for item in rotation_list:
        #     print(item.text)

        # 生成 suffix array 和 BWT
        for item in rotation_list:
            suffix_array.append(item.pos)
            bwt.append(item.text[-1:])

        # 生成 reference 的 Occ
        for item in rotation_list:
            for symbol in alphabet:
                if len(Occ[symbol]) == 0:
                    prev = 0
                else:
                    prev = Occ[symbol][-1]

                if item.text[-1:] == symbol:
                    Occ[symbol].append(prev + 1)
                else:
                    Occ[symbol].append(prev)

        self.data = dict()
        self.data['S'] = suffix_array
        self.data['B'] = bwt
        self.data['C'] = C
        self.data['O'] = Occ
        self.text = reference
        self.text_num = len(reference_array)
    
    # 生成反向互补序列
    def add_reverse_reference(self, reference):

        reference_list = list()

        for item in reference:
            reference_list.append(item)
            item_reverse = ''.join([self.reverse_dict[char] for char in list(item)][::-1])
            reference_list.append(item_reverse)

        return reference_list

    def Occ(self, char, index):
        if index < 0:
            return 0
        else:
            return self.data['O'][char][index]
