import time
import argparse
import sys

# 通过命令行传递参数
parser = argparse.ArgumentParser()
parser.add_argument('--difference_threshold', '-d', type=int, required=False)  # 自定义差异数值上限
parser.add_argument('--reference_file', '-r', type=str, required=True)         # 参考序列文件
parser.add_argument('--query_file', '-q', type=str, required=True)             # 待比对序列文件

args = parser.parse_args()

difference_threshold = args.difference_threshold if args.difference_threshold is not None else 2 # 默认差异数值上限2
reference_file = args.reference_file
query_file = args.query_file

# 通过文件读取参考序列和待比对序列
reference = ''
query = ''


# 文件读取
with open(reference_file) as file_obj:
    for content in file_obj:
        reference = reference + content.replace(' ', '').replace('\n', '').replace('\r', '')


with open(query_file) as file_obj:
    for content in file_obj:
        query = query + content.replace(' ', '').replace('\n', '').replace('\r', '')


class BWA:
    # 初始化函数生成所有和reference相关的数据结构
    def __init__(self, reference):
        """
        :param reference: 参考基因组序列
        """
        reference = reference.lower()

        rotation_list, rotation_list_reverse, suffix_array, suffix_array_reverse, bwt = [list() for i in range(5)]
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

        # 生成 suffix_array_reverse
        for item in rotation_list_reverse:
            suffix_array_reverse.append(item.pos)

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
        self.C = C
        self.Occ = Occ
        self.Occ_reverse = Occ_reverse
        self.n = len(reference)  # 参考基因组序列长度
        self.D = list()  # lower bound on the number of differences
        self.alphabet = alphabet

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


        # 两个递归的中止条件
        # stop1：允许的差异数过少，无法满足匹配要求，即 z < D(i)
        if z < self.get_D(i):
            return set()  # 返回空集
        # stop2：整个短序列都被搜索完成，返回 SA 的索引的上下界
        if i < 0:
            for m in range(k, l + 1):
                tempset.add(m)
            return tempset

        result = set()

        result = result.union(self.InexRecur(read, i - 1, z - 1, k, l))  # 插入

        for symbol in self.alphabet:
            newK = self.C[symbol] + self.OCC(symbol, k-1) + 1
            newL = self.C[symbol] + self.OCC(symbol, l)
            if newK <= newL:
                # 单个字符找到
                result = result.union(self.InexRecur(read, i, z - 1, newK, newL))  # 删除
                if symbol == read[i]:
                    # 该字符成功比对, 继续比对但不减少 z
                    result = result.union(self.InexRecur(read, i - 1, z, newK, newL))
                else:
                    # 该字符匹配失败，继续比对但减少 z，意味着一次不匹配
                    result = result.union(self.InexRecur(query, i - 1, z - 1, newK, newL))

        return result

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
    def match(self, read, z):
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


class Suffix:
    def __init__(self, text, position):
        self.text = text
        self.pos = position


def textKey(a): return a.text


if __name__ == "__main__":

    # 参考序列和待比对序列不存在/长度为0
    if len(reference) == 0 or len(query) == 0:
        print("Reference or Query DON't EXIST")
        sys.exit(0)

    time_start = time.time()  # 开始计时

    data = BWA(reference)
    print("\nReference: \"%s\"" % reference)

    print("Searching for \"%s\" with max difference threshold of %d..." % (query, difference_threshold))

    matches = data.match(query, difference_threshold)

    print("\n%d match(es) at position(s): %s \n\n" % (len(matches), matches))

    time_end = time.time()  # 结束计时
    time_c = time_end - time_start  # 运行所花时间

    print('耗时:', time_c, 's\n')

    # 将比对结果输出到 .sai.txt 文件

    f = open('ref_query.sai.txt', 'w')
    for match in matches:
        f.write(str(match)+',')
    f.write(str(difference_threshold))
    f.close()










