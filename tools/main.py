import bwt
import algorithm

reference =["ACG", "GTG"] # 不包含终止符 $


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

    print(algorithm.get_seq(2, fmd_index.data['B'], fmd_index.data['O'], fmd_index.data['C']))
    print(algorithm.get_seq(3, fmd_index.data['B'], fmd_index.data['O'], fmd_index.data['C']))