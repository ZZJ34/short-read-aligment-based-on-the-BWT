"""
Burrows-Wheeler Alignment - String Search
--------------------------------------------
This is a very simply quick and dirty implementation of a Burrows-Wheeler Aligner for indexing and sequence alignment.
It does NOT implement any of the heuristic methods included in the actual BWA algorithm, so implements the basic algorithm which returns 100% accurate results
It uses the Burrows-Wheeler transform (BWT), Suffix Array (SA), and 2 other auxillary datastructures C and Occ (sometimes called O)
It uses the D array to prune the tree search and speed up the inexact search algorithm.
The search is case insensitive.
Differences between this code and the real BWA algorithm
--------------------------------------------------------
- This is NOT IN ANY WAY optimised. This has been coded for ease of understanding, ease of reading and with the goal of better understaning the basic BWA algorithm and its datastructures
- This code parses the suffix tree using a recursive depth-first search, while the real BWA uses breadth-first search (and a heap datastructure), and uses this to prioritise the best partial alignments first
- If BWA finds a result with a difference score of z, it only further considers matches with difference z+1, speeding up the process by ignoring worse results and pruning the search space more aggresively
- BWA sets a maximum allowed difference in the first few tens of bases (called the seed) resulting in 2.5x improvement, and very little loss in accuracy. This is not effective on shorter reads (<50bp)
- BWA reduces the required operating memory by storing small fractions of the Occ and SA arrays, and calculating the rest on the fly. This implementation does not do this.
"""

"""
并非自己原创，参考Jwomers的代码进行修改。
源代码：https://github.com/Jwomers/burrows_wheeler_alignment
原文章：Fast and accurate short read alignment with Burrows-Wheeler transform
"""

import time
import math
import traceback

# environment variables
debug = True
only_InexRecur = True  # 强制使用 InexRecur 递归过程实现精确/非精确匹配
use_find_match = True  # 使用 find_match 或者 find_match_2
show_data_structures = False
use_lower_bound_tree_pruning = True  # set this to false (in conjunction with debug=True) to see the full search through the suffix trie
show_data_compress = False  # 是否展示数据压缩的过程
"""
根据文章所诉，使用 lower_bound 可以有效减少搜索空间
"""
# search parameters
indels_allowed = True# turn off for mismatches only, no insertion or deletions allowed
difference_threshold = 1
insertion_penalty = 1
deletion_penalty = 1
mismatch_penalty = 1
# compress parameters
block_number = 2  # 分成的块数，按照文章所述的内容，仅仅分成两块
use_middle_as_head_number = False  # 是否使用分块中间的条目作为标头，否就用第一个条目

# reference and query strings
reference ="""CCTGAG"""
query = "CTG"

"""
A Burrows-Wheeler Alignment class
---------------------------------
包含了四种该算法涉及的核心数据结构，SA（Suffix array/后缀数组），BWT string（BWT字符串），C，Occ（也称为O）
---------------------------------
关于数据结构C和Occ请参考原文章的2.4节
"""
class BWA:
    # 初始化函数生成所有和reference相关的数据结构
    def __init__(self, reference):
        """
        :param reference: 参考基因组序列
        """
        reference = reference.lower()

        rotation_list, rotation_list_reverse, suffix_array, suffix_array_reverse, bwt, bwt_reverse = [list() for i in range(6)]
        C, Occ, Occ_reverse = [dict() for i in range(3)]
        alphabet = set()  # 创建一个无序不重复元素集
        reverse_reference = reference[::-1]  # 参考基因序列的反向序列

        # 创建字符表
        for char in reference:
            alphabet.add(char)  # 此时 alphabet 为'g', 'a', 't', 'c'的集合


        # 初始化C,Occ,Occ_reverse
        for char in alphabet:
            C[char] = 0
            Occ[char] = list()  # in Occ, each character has an associated list of integer values (for each index along the reference)
            Occ_reverse[char] = list()

        # 将结束字符添加到参考基因组字符串
        reference = "%s$" % reference
        reverse_reference = "%s$" % reverse_reference

        # 生成基因组参考序列的循环移位序列，并存储于 rotation_list
        for i in range(len(reference)):
            new_rotation = "%s%s" % (reference[i:], reference[0:i])
            struct = Suffix(new_rotation, i)

            rotation_list.append(struct)
        # 生成逆序基因组参考序列的循环移位序列，并存储于 rotation_list_reverse
        for i in range(len(reverse_reference)):
            new_rotation_reverse = "%s%s" % (reverse_reference[i:], reverse_reference[0:i])
            struct_rev = Suffix(new_rotation_reverse, i)

            rotation_list_reverse.append(struct_rev)

        # 生成C
        # C(a) = the number of characters 'a' in the Reference that are lexographically smaller than 'a'
        for symbol in alphabet:
            for i in range(len(reference)-1):
                if reference[i] < symbol:
                    C[symbol] = C[symbol] + 1

        # 将循环移位序列根据 text 的字典序排序
        rotation_list.sort(key=textKey)
        rotation_list_reverse.sort(key=textKey)

        # 生成 suffix array 和 BWT
        for item in rotation_list:
            suffix_array.append(item.pos)
            bwt.append(item.text[-1:])

        # 生成 suffix_array_reverse 和 BWT_reverse
        for item in rotation_list_reverse:
            suffix_array_reverse.append(item.pos)
            bwt_reverse.append(item.text[-1:])

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

        # 生成 reverse_reference 的 Occ
        for item in rotation_list_reverse:
            for symbol in alphabet:
                if len(Occ_reverse[symbol]) == 0:
                    prev = 0
                else:
                    prev = Occ_reverse[symbol][-1]

                if item.text[-1:] == symbol:
                    Occ_reverse[symbol].append(prev + 1)
                else:
                    Occ_reverse[symbol].append(prev)

        # save all the useful datastructures as class variables for easy future access
        self.SA = suffix_array
        self.SA_reverse = suffix_array_reverse
        self.BWT = bwt
        self.BWT_reverse = bwt_reverse
        self.C = C
        self.Occ = Occ
        self.Occ_reverse = Occ_reverse
        self.n = len(reference)  # 参考基因组序列长度
        self.D = list()  # lower bound on the number of differences
        self.alphabet = alphabet
        self.InexRecur_parameter_list = []  # 递归调用 InexRecur 过程中的参数列表

    # 工具：计算 Occ(a,i) 和 Occ_reverse(a,i) 的值，通过 reverse 判断 Occ 还是 Occ_reverse
    # 注意：Occ(a,-1) = 0
    def OCC(self, char, index, reverse=False):
        if index < 0:
            return 0
        else:
            if reverse:
                return self.Occ_reverse[char][index]
            else:
                return self.Occ[char][index]

    # 工具：计算短序列的 D，用于修建遍历树并提升搜寻的速度
    def calculated_D(self, read):
        read = read.lower()

        tempset = set()

        k = 0
        l = self.n - 1
        z = 0

        for i in range(len(read)):
            k = self.C[read[i]] + self.OCC(read[i], k - 1, reverse=True) + 1
            l = self.C[read[i]] + self.OCC(read[i], l, reverse=True)
            if k > l:
                k = 0
                l = self.n - 1
                z = z + 1
            self.D.append(z)

        for m in range(k, l + 1):
            tempset.add(m)
        return tempset

    # 工具：获取 D 中的数值
    def get_D(self, index):
        if index < 0:
            return 0
        else:
            return self.D[index]

    # 工具：非精确搜索的递归函数
    def InexRecur(self, read, i, z, k, l):
        tempset = set()

        self.InexRecur_parameter_list.append([i, z, k, l])

        # 两个递归的中止条件
        # stop1：允许的差异数过少，无法满足匹配要求，即 z < D(i)
        if (use_lower_bound_tree_pruning and z < self.get_D(i)) or (not use_lower_bound_tree_pruning and z < 0):
            # 不使用 D 时，就默认 D(i) = 0
            if debug:
                print("too many differences, terminating path\n")
            return set()  # 返回空集
        # stop2：整个短序列都被搜索完成，返回 SA 的索引的上下界
        if i < 0:
            if debug:
                print("query string finished, terminating path, success! k=%d, l=%d\n" % (k, l))
            for m in range(k, l + 1):
                tempset.add(m)
            return tempset

        result = set()
        if indels_allowed:
            # 允许插入和删除
            result = result.union(self.InexRecur(read, i - 1, z - insertion_penalty, k, l))  # 插入

        for symbol in self.alphabet:
            newK = self.C[symbol] + self.OCC(symbol, k-1) + 1
            newL = self.C[symbol] + self.OCC(symbol, l)
            if newK <= newL:
                # 单个字符找到
                if indels_allowed:
                    # 允许插入和删除
                    result = result.union(self.InexRecur(read, i, z - deletion_penalty, newK, newL))  # 删除
                if debug:
                    print("char '%s found' with k=%d, l=%d. z = %d: parent k=%d, l=%d" % (symbol, newK, newL, z, k, l))
                if symbol == read[i]:
                    # 该字符成功比对, 继续比对但不减少 z
                    result = result.union(self.InexRecur(read, i - 1, z, newK, newL))
                else:
                    # 该字符匹配失败，继续比对但减少 z，意味着一次不匹配
                    result = result.union(self.InexRecur(query, i - 1, z - mismatch_penalty, newK, newL))

        return result

    def match(self, read, num_differences):
        if use_find_match:
            if debug:
                print("Exact Matching + Inexact Matching\n")
            return self.find_match(read, num_differences)
        else:
            if debug:
                print("Exact Matching/Inexact Matching\n")
            return self.find_match_2(read, num_differences)

    # 查找短序列在参考基因中的位置：精确匹配 + 非精确匹配
    def find_match(self, read, num_differences):
        """
        :param read: 查找的短序列
        :param num_differences: 最大允许的差异数（SNP，插入，删除）
        :return: 匹配的位置
        """
        """
        此处要说明的：
        inexact_match 并非只能完成非精确匹配 
        当 num_differences 为 0 也能完成精确匹配，但效率略低
        """
        if num_differences == 0 and not only_InexRecur:
            # 精确匹配（backward search/向后搜索）
            return self.exact_match(read)
        else:
            # 非精确匹配（bounded traversal/边界遍历）
            return self.inexact_match(read, num_differences)

    # 查找短序列在参考基因中的位置：精确匹配/非精确匹配
    def find_match_2(self, read, num_differences):
        return self.all_match(read, num_differences)

    # 精确匹配算法
    def exact_match(self, read):
        """
        :param read: 查找的短序列
        :return: 匹配的位置
        """
        read = read.lower()

        # 对于空串，i = 0，j = n-1
        # 从空串开始搜寻
        k = 0
        l = self.n - 1

        for x in range(len(read)):
            newRead = read[-1-x]
            newK = self.C[newRead] + self.OCC(newRead, k - 1) + 1
            newL = self.C[newRead] + self.OCC(newRead, l)
            k = newK
            l = newL
        matches = self.SA[k:l + 1]  # 这里是因为 python 的数组截取是左闭右开，k 是下边界，l 是上边界
        return matches

    # 非精确匹配算法
    def inexact_match(self, read, z):
        """
        :param read: 查询的短序列
        :param z: 允许的最大不匹配数
        :return: 匹配结果
        """
        read = read.lower()
        self.calculated_D(read)
        SA_index = self.InexRecur(read, len(read)-1, z, 0, self.n - 1)
        return [self.SA[x] for x in SA_index]

    # Exact Matching Using Reversed Reference Genome data
    # 同时实现精确/非精确匹配
    # 该算法的实现参考：Fast and accurate short read alignment with Burrows-Wheeler transform 中的算法2
    def all_match(self, read, z):
        read = read.lower()

        # 计算
        I_reverse = self.calculated_D(read)

        # 判断：如果 D 所有的元素都是则可以实现精确匹配
        sum_D = 0
        for item in self.D:
            sum_D = sum_D + item

        if sum_D == 0:
            # 精确匹配
            """
            #这里要说明一下：
            --------------------
            文献里这里没有减1，但是算法实现的时候需要减1，才能和另一种方法的结果吻合
            个人猜测，可能是对参考序列长度的定义有不同，即是否包含最后的$
            """
            return [self.n - self.SA_reverse[x] - len(read) - 1 for x in I_reverse]
        else:
            # 非精确匹配
            SA_index = self.InexRecur(read, len(read) - 1, z, 0, self.n - 1)
            return [self.SA[x] for x in SA_index]


"""
原文章 Hardware-Acceleration of Short-Read Alignment Based on the Burrows-Wheeler Transform
------------------------------------------------------------------------------------------
根据文章中所诉的方法，对于 Occ 辅助数据结构进行压缩已经解码
"""
class OCC_COMPRESS:
    def __init__(self, Occ, bwt_str):
        # 确定块的长度
        block_len = math.ceil(len(bwt_str) / block_number)
        # 确定标头的位置
        if use_middle_as_head_number:
            head_index = math.ceil(block_len/2) - 1
        else:
            head_index = 0

        # 压缩后的数据初始化
        compressed_data = {}
        for i in range(block_number):
            compressed_data[i] = []

        # 符号排序
        alphabet = list(Occ.keys())
        alphabet.sort()

        # 不包含最后一个 code
        for i in range(block_number):
            for item in alphabet:
                compressed_data[i].append(Occ[item][i * block_len + head_index])
            compressed_data[i].extend(bwt_str[1 + i * block_len:(i + 1) * block_len])

        self.Occ = Occ                          # 原始数据
        self.block_len = block_len              # 块长
        self.head_index = head_index            # 标头位置索引
        self.alphabet = alphabet                # 标头符号的顺序
        self.compressed_data = compressed_data  # 压缩的数据

    # 根据压缩数据进行解码
    def decode(self, char, index):
        # 判断 index 的索引范围
        if index > (len(self.Occ[char]) - 1):
            raise Exception("list index out of range")
        if index < 0:
            return 0

        # 判断标头的位置
        if use_middle_as_head_number:
            code_index = math.floor(index / self.block_len)
            code_item = self.compressed_data[code_index]
            if index - code_index * self.block_len < self.head_index:  # 位于标头之前
                count_str = code_item[4:][index - code_index * self.block_len:self.head_index]
                char_index = self.alphabet.index(char)

                # head 的数值
                count = code_item[char_index]

                # 统计 body 中的 char 数量
                for i in count_str:
                    if i == char: count = count - 1

                return count
            if index - code_index * self.block_len == self.head_index:  # 位于标头
                return code_item[self.alphabet.index(char)]
            if index - code_index * self.block_len > self.head_index:  # 位于标头之后

                count_str = code_item[4:][self.head_index:index - code_index * self.block_len]
                char_index = self.alphabet.index(char)

                # head 的数值
                count = code_item[char_index]

                # 统计 body 中的 char 数量
                for i in count_str:
                    if i == char: count = count + 1

                return count

        else:
            code_index = math.floor(index / self.block_len)
            code_item = self.compressed_data[code_index]

            count_str = code_item[4:][:index - code_index * self.block_len]

            char_index = self.alphabet.index(char)

            # head 的数值
            count = code_item[char_index]

            # 统计 body 中的 char 数量
            for i in count_str:
                if i == char: count = count + 1

            return count
"""
A simple class with 2 variables, used for sorting and calculating the Suffix Array and BWT array
Each instance holds the position of the suffix and the suffix (text) itself
"""
class Suffix:
    def __init__(self, text, position):
        self.text = text
        self.pos = position

"""
this is used to sort the Suffix objects, according to their text key
"""
def textKey( a ): return a.text


if __name__ == "__main__":

    time_start = time.time()  # 开始计时

    data = BWA(reference)
    print("\n\nReference: \"%s\"" % reference)
    print("length of reference: %i\n" % len(reference))
    if show_data_structures:
        # printing out the datastructues for manual inspection
        print("\nSA		BWT-Str")
        print("--		-------")
        for i in range(len(data.SA)):
            print("%s		   %s" % (data.SA[i], data.BWT[i]))
        print("\n")
        print("C(a) = the number of characters 'a' in the Reference that are lexographically smaller than 'a'")
        print(data.C)
        print("\n")
        print("Occ(a,i) is the number of occurances of the character 'a' in BWT[0,i]")
        for str in data.alphabet:
            print(str, end="  ")
            print(data.Occ[str])
        print("\n")
        print("Occ_reverse(a,i) is the number of occurances of the character 'a' in BWT_reverse[0,i]")
        for str in data.alphabet:
            print(str, end="  ")
            print(data.Occ_reverse[str])
        print("\n")

    if indels_allowed:
        extra = " and insertions and deletions allowed"
    else:
        extra = " with no insertions or deletions allowed"
    print("Searching for \"%s\" with max difference threshold of %d%s...\n" % (query, difference_threshold, extra))

    matches = data.match(query, difference_threshold)

    if show_data_structures:
        if use_lower_bound_tree_pruning:
            print("\nD array:", end=" ")
            print(data.D)
        else:
            print("\nThere is no calculation about the D array for a query")

    print("\n%d match(es) at position(s): %s \n\n" % (len(matches), matches))

    time_end = time.time()  # 结束计时

    time_c = time_end - time_start  # 运行所花时间
    print('耗时:', time_c, 's\n')

    if show_data_structures:
        print("展示 InexRecur 递归调用过程中的参数:")
        print("index  i   z   k   l")
        print("-----  -   -   -   -")
        index = 1
        for item in data.InexRecur_parameter_list:
            print("%3d  %3d %3d %3d %3d" % (index, item[0], item[1], item[2], item[3]))
            index = index + 1
        print("\n")

    if show_data_compress:
        # 展示对 Occ 的压缩和解码
        occ_compress = OCC_COMPRESS(data.Occ, data.BWT)

        print("\n针对 Occ 压缩解码的实验比对")
        print("解码结果：Occ(%s,%i)=%i" % ("t", 16, occ_compress.decode("t", 16)))
        print("原码结果：Occ(%s,%i)=%i\n" % ("t", 16, data.OCC("t", 16)))

        # 展示 Occ 压缩数据
        if show_data_structures:
            print('压缩后的Occ')
            for i in range(block_number):
                print("code %i" % i, end=" ")
                print(occ_compress.compressed_data[i])

        occ_reverse_compress = OCC_COMPRESS(data.Occ_reverse, data.BWT_reverse)

        print("\n针对 Occ_reverse 压缩解码的实验比对")
        print("解码结果：Occ(%s,%i)=%i" % ("t", 3, occ_reverse_compress.decode("t", 3)))
        print("原码结果：Occ(%s,%i)=%i\n" % ("t", 3, data.OCC("t", 3, reverse=True)))

        # 展示 Occ_reverse 压缩数据
        if show_data_structures:
            print('压缩后的Occ_reverse')
            for i in range(block_number):
                print("code %i" % i, end=" ")
                print(occ_reverse_compress.compressed_data[i])




















