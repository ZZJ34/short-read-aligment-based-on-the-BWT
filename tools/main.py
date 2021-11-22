import bwt

reference =["ACG", "GTG", "ATTG"] # 不包含终止符 $


if __name__ == "__main__":
    print('hello, my bwa tools!')

    fm_index = bwt.BWA_FM_index(reference)


    print('********************** RESULT *************************')

    print('B :', fm_index.data['B'])
    print('S :', fm_index.data['S'])
    print('C :', fm_index.data['C'])
    print('O :', fm_index.data['O'])

    print(fm_index.text)