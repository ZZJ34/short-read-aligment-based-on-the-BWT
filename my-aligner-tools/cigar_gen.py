
import numpy as np
import argparse

# 通过命令行传递参数
parser = argparse.ArgumentParser()
parser.add_argument('--sai_file', '-s', type=str, required=True)               # 比对位置
parser.add_argument('--reference_file', '-r', type=str, required=True)         # 参考序列文件
parser.add_argument('--query_file', '-q', type=str, required=True)             # 待比对序列文件

args = parser.parse_args()

reference_file = args.reference_file
query_file = args.query_file
sai_file = args.sai_file

# 通过文件读取参考序列和待比对序列
reference_original = ''
query_original = ''
matches = list()
difference_threshold = 0

# 文件读取
with open(reference_file) as file_obj:
    for content in file_obj:
        reference_original = reference_original + content.replace(' ', '').replace('\n', '').replace('\r', '')

with open(query_file) as file_obj:
    for content in file_obj:
        query_original = query_original + content.replace(' ', '').replace('\n', '').replace('\r', '')

with open(sai_file) as file_obj:
    content = file_obj.read().split(',')
    for i in range(len(content)-1):
       matches.append(int(content[i]))
    difference_threshold = int(content[len(content)-1])

# 环境变量
TYPE_MATCH = 0      # 匹配
TYPE_MISMATCH = 1   # SNP
TYPE_INSERTION = 2  # 相对参考序列插入，待比对序列比参考序列多一个
TYPE_DELETION = 3   # 相对参考序列删除，待比对序列比参考序列少一个
# 不同匹配类型的罚分
# 不同的罚分值可能会得到不同的 CIGAR
PENALTY_MATCH = 5
PENALTY_MISMATCH = -2
PENALTY_INSERTION = -4
PENALTY_DELETION = -4

# M:匹配 X：SNP I：插入 D：删除
kTypeToSymbol = ['M', 'X', 'I', 'D']


class CIGAR:
    def __init__(self, reference, query, postion):
        self.reference = reference.upper()
        self.query = query.upper()

        self.cigar = ""
        self.edit_distance = 0
        self.alignment = ""
        self.max_score = 0
        self.num_insertions = 0
        self.num_deletions = 0
        self.position = postion


    def generateCigar(self):
        # 参考序列和待比对序列不存在/长度为0
        if self.reference is None or self.query is None or len(self.reference) == 0 or len(self.query) == 0:
               return -1

        # DP matrix and traceback, and initialize all values to 0
        dp_matrix = np.zeros((len(self.query)+1, len(self.reference)+1))
        dp_traceback = np.zeros((len(self.query)+1, len(self.reference)+1), dtype=np.int)

        up, left, diagonal = [0, 0, 0]

        for i in range(len(self.query) + 1):
            dp_matrix[i][0] = i * PENALTY_INSERTION
        for i in range(len(self.reference) + 1):
            dp_matrix[0][i] = i * PENALTY_DELETION

        # 记录最大评分数值以及索引
        max_score_index = 0

        # 填充评分矩阵
        for i in range(1, len(self.query) + 1):
            for j in range(1, len(self.reference) + 1):
                up = dp_matrix[i - 1][j] + PENALTY_INSERTION
                left = dp_matrix[i][j - 1] + PENALTY_DELETION
                diagonal = dp_matrix[i - 1][j - 1] + (PENALTY_MATCH if self.query[i-1] == self.reference[j-1] else PENALTY_MISMATCH)

                # 判断三者当中的最大值
                dp_matrix[i][j] = diagonal
                dp_traceback[i][j] = 0 if self.query[i - 1] == self.reference[j - 1] else 1
                if up > dp_matrix[i][j]:
                    dp_matrix[i][j] = up
                    dp_traceback[i][j] = 2
                if left > dp_matrix[i][j]:
                    dp_matrix[i][j] = left
                    dp_traceback[i][j] = 3

                # 判断最后一行的最大值
                if i == len(self.query) and dp_matrix[len(self.query)][j] > dp_matrix[len(self.query)][max_score_index]:
                    max_score_index = j
                    self.max_score = dp_matrix[len(self.query)][j]

        # 回溯整个矩阵
        current_traceback_row = len(self.query)
        current_traceback_column = max_score_index

        # 找到 CIGAR 字符串
        iterations = 0
        type = 0
        last_type = 0
        num_same_types = 0

        while current_traceback_row > 0 and current_traceback_column > 0:
            type = dp_traceback[current_traceback_row][current_traceback_column]

            # 记录回溯 DP_matrix 向同一个方向移动的次数（上，左以及对角线）
            if iterations == 0:
                # 第一次开始回溯
                last_type = type
                num_same_types = 1
            else:
                # 比对结果没有变
                if type == last_type:
                    num_same_types = num_same_types + 1
                # 比对结果发生改变
                else:
                    self.cigar = str(num_same_types) + kTypeToSymbol[last_type] + self.cigar
                    last_type = type
                    num_same_types = 1

            # 统计 edit_distance
            if not type == TYPE_MATCH:
                self.edit_distance = self.edit_distance + 1

            # 统计插入和删除数量
            if type == TYPE_DELETION:
                self.num_deletions = self.num_deletions + 1
            if type == TYPE_INSERTION:
                self.num_insertions = self.num_insertions + 1

            # 给出 alignment
            if type == TYPE_DELETION:
                self.alignment = '-' + self.alignment                     # 相对参考序列缺失（少一个）
            elif type == TYPE_MATCH:
                self.alignment = self.query[current_traceback_row - 1] + self.alignment    # 相对参考序列比对成功
            elif type == TYPE_INSERTION:
                self.alignment = '+' + self.alignment                      # 相对参考序列插入（多一个）
            else:
                self.alignment = 'x' + self.alignment                     # 相对查考序列 SNP

            current_traceback_row = current_traceback_row - (1 if not type == TYPE_DELETION else 0)

            current_traceback_column = current_traceback_column - (1 if not type == TYPE_INSERTION else 0)

            iterations = iterations + 1

        # 将最后剩余的 CIGAR 统计加入
        self.cigar = str(num_same_types) + kTypeToSymbol[last_type] + self.cigar


if __name__ == "__main__":

    cigar_list = list()

    for match in matches:
        i = 0
        cigar_temp = CIGAR(reference_original[match - i:], query_original, match)
        cigar_temp.generateCigar()
        # while len(cigar_temp.alignment) - cigar_temp.num_deletions < len(query_original):
        #     i = i + 1
        #     cigar_temp = CIGAR(reference_original[match - i:], query_original, match)
        #     cigar_temp.generateCigar()

        cigar_list.append(cigar_temp)

    for item in cigar_list:
        print(item.edit_distance, end=' ')
        print(item.max_score, end=' ')
        print(item.position, end=' ')
        print(item.cigar, end=' ')
        print(item.alignment, end=' ')
        print('\n')



