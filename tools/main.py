import bwt
import algorithm
import copy

# $, A, C, G, T, N/X 从 0 开始编码
# 前向扩展 => W -> Wa
# 后向扩展 => W -> aW
# 对于空串，k = l = 0, s = len(T)-1


reference =["ACGTX"] # 不包含终止符 $


if __name__ == "__main__":
    print('hello, my bwa tools!')

    fm_index = bwt.BWA_FM_index(reference)

    fmd_index = bwt.BWA_FMD_index(reference)


    print('********************** RESULT *************************')

    # print('B :', fm_index.data['B'])
    # print('S :', fm_index.data['S'])
    # print('C :', fm_index.data['C'])
    # print('O :', fm_index.data['O'])

    # print(fm_index.text)

    print('B :', fmd_index.data['B'])
    print('S :', fmd_index.data['S'])
    print('C :', fmd_index.data['C'])
    print('O :', fmd_index.data['O'])

    print(fmd_index.text)

    print('********************** Sequence retrieval *************************')

    # print(algorithm.get_all_seq(len(reference)*2, fmd_index))
    
    # print(algorithm.backward_extension([0, 0, len(fmd_index.text)-1], 'A', fmd_index))
    # print(algorithm.forward_extension([0, 0, len(fmd_index.text)-1], 'A', fmd_index))

    print(algorithm.super_MEM1("ACG", 1, fmd_index))

    