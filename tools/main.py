import bwt
import data
import algorithm

# $, A, C, G, T, N/X 从 0 开始编码
# 前向扩展 => W -> Wa
# 后向扩展 => W -> aW
# 对于空串，k = l = 0, s = len(T)-1


# reference =["ACCTTGA"] # 不包含终止符 $


if __name__ == "__main__":
    print('hello, my bwa tools!')


    
    print('********************** DATA LOAD *************************')

    reference = data.load_file('d:\\short-read-aligment-based-on-the-BWT\\tools\\my_ref.fa')
    
    # fmd_index = bwt.BWA_FMD_index(reference)

    fmd_index = bwt.BWA_FMD_index_noend(reference)

    print('************************ RESULT **************************')

    # print('B :', fm_index.data['B'])
    # print('S :', fm_index.data['S'])
    # print('C :', fm_index.data['C'])
    # print('O :', fm_index.data['O'])

    # print(fm_index.text)

    # print('B   :', fmd_index.data['B'])
    # print('B_n :', fmd_index.data['B_n'])
    # print('S   :', fmd_index.data['S'])
    # print('C   :', fmd_index.data['C'])
    # print('O   :', fmd_index.data['O'])
    # print('O_n :', fmd_index.data['O_n'])

    # print(fmd_index.text)

    # print(len(fmd_index.data['B_n']))

    print('********************** ALGORITHM *************************')

    # print(algorithm.get_all_seq(len(reference)*2, fmd_index))
    
    # print(algorithm.backward_extension([0, 0, len(fmd_index.text)-1], 'A', fmd_index))
    # print(algorithm.forward_extension([0, 0, len(fmd_index.text)-1], 'A', fmd_index))

    # print(algorithm.super_MEM1("ACTTG", 0, fmd_index))

    data.data_store(fmd_index)


    data.data_store_encode(fmd_index)



    