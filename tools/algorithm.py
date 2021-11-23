
symbols = ['$', 'A', 'C', 'G', 'T', 'X']

reverse_dict={
    "A" : "T",
    "T" : "A",
    "C" : "G",
    "G" : "C",
    "X" : "X"
}

# Algorithm 1: Sequence retrieval
# Input: Sequence index i≥0; B, O and C defined in the text
# Output: Sequence P and k, the rank of P
def get_seq(index, all_data):
    k = index
    p = ''

    while True:
        a = all_data.data['B'][k]
        k = all_data.data['C'][a] + all_data.Occ(a, k) - 1
        
        if(a == '$'): break

        p = a + p

    return (p, k)

def get_all_seq(reference_num, all_data):
    p_k_list = []
    # 既然 get_seq 的 index 不准确，那就全部获取
    for i in range(reference_num):
        (p,k) = get_seq(i, all_data)
        p_k_list.append((p,k))
    return p_k_list

# Algorithm 2: Backward extension
# Input: Bi-interval [k,l,s] of string W and a symbol a
# Output: Bi-interval of string aW
def backward_extension(interval, symbol, all_data):
    [k, l, s] = interval
    k_list, s_list, l_list= [list() for i in range(3)]  
    for item in symbols:
        k_list.append(all_data.data['C'][item] + all_data.Occ(item, k-1))
        s_list.append(all_data.Occ(item, k+s-1) - all_data.Occ(item, k-1))

    l_list = [0 for i in range(len(symbols))]

    l_list[0] = l
    l_list[4] = l_list[0] + s_list[0]
    l_list[3] = l_list[4] + s_list[4]
    l_list[2] = l_list[3] + s_list[3]
    l_list[1] = l_list[2] + s_list[2]
    l_list[5] = l_list[1] + s_list[1]

    return [k_list[symbols.index(symbol)], l_list[symbols.index(symbol)], s_list[symbols.index(symbol)]]
    
# Algorithm 3: Forward extension
# Input: Bi-interval [k,l,s] of string W and a symbol a
# Output: Bi-interval of string Wa
def forward_extension(interval, symbol, all_data):
    return backward_extension(interval, reverse_dict[symbol], all_data)


        


