
symbols = ['A', 'C', 'G', 'T']

symbols_encode = {
    'A' : '00', 
    'C' : '01', 
    'G' : '10', 
    'T' : '11'
}

def load_file(file_dir):
    line_list = list()
    with open(file_dir,'r') as f:
        lines=f.readlines()
        for line in lines:
            if line.startswith('>'):
                line_list.append(list())
            else:
                line_list[-1].append(line.split('\n')[0])

    return [''.join(reference).replace('N', 'X') for reference in line_list]

def data_store(all_data, interval=128):
    # 明码存储
    data = []  
    i = 1
    while True:
        if interval*i < len(all_data.data['B_n']):
            Occ_str = ''
            for symbol in symbols:
                Occ_str = Occ_str + "Occ(%s, %d)=%d " % (symbol, interval*i-1, all_data.Occ_n(symbol,interval*i-1))
        
            BWT_str = ''.join(all_data.data['B_n'][(i-1)*interval:i*interval])
            data.append(Occ_str + BWT_str)
        else:
            Occ_str = ''
            for symbol in symbols:
                Occ_str = Occ_str + "Occ(%s, %d)=%d " % (symbol, len(all_data.data['B_n'])-1, all_data.Occ_n(symbol, len(all_data.data['B_n'])-1))
        
            BWT_str = ''.join(all_data.data['B_n'][(i-1)*interval:])
            data.append(Occ_str + BWT_str)
            break

        i = i+1

    f = open('./bwt.txt','w')
    f.write(str(all_data.primary)+'\n')
    for data_line in data:
        f.write(data_line+'\n')
    f.close()

    return data


def data_store_encode(all_data, interval=128):
    # 编码存储
    data = []  
    i = 1
    while True:
        if interval*i < len(all_data.data['B_n']):
            Occ_str = ''
            for symbol in symbols:
                Occ_str = Occ_str + '{:064b}'.format(all_data.Occ_n(symbol,interval*i-1))
            
            BWT_str = ''.join([ symbols_encode[item] for item in all_data.data['B_n'][(i-1)*interval:i*interval]])
            data.append(Occ_str + BWT_str)
        else:
            Occ_str = ''
            for symbol in symbols:
                Occ_str = Occ_str + '{:064b}'.format(all_data.Occ_n(symbol, len(all_data.data['B_n'])-1))
            
            BWT_str = ''.join([ symbols_encode[item] for item in all_data.data['B_n'][(i-1)*interval:]])
            data.append(Occ_str + BWT_str)
            break

        i = i+1
    


    f = open('./bwt_encode.txt','w')
    f.write('{:064b}'.format(all_data.primary)+'\n')
    for data_line in data:
        f.write(data_line+'\n')
    f.close()