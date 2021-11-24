
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


def data_store(name, all_data):
    # 明码存储
    print(name)

def data_store_encode(name, all_data):
    # 编码存储
    print(name)