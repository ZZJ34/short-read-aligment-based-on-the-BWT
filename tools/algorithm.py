import copy

# 按照字典序排序


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

# Algorithm 5: Finding SMEMs
# Input: String P and start position i0; 特殊定义 P[−1]=0, P[−1]=$
# Output: Set of bi-intervals of SMEMs overlapping i0
def super_MEM1(P, start_index, all_data):

    # 数据准备
    curr_list, prev_list, match_list = [list() for i in range(3)]

    # 初始数据
    [k, l, s] = [
        all_data.data['C'][P[start_index]], 
        all_data.data['C'][reverse_dict[P[start_index]]],
        all_data.data['C'][symbols[symbols.index(P[start_index])+1]] - all_data.data['C'][P[start_index]]
    ]

    # 前向搜索
    i = start_index + 1
    while i <= len(P):
        if i == len(P):
            curr_list.append([k, l, s]) 
        else:
            [k_temp, l_temp, s_temp] = forward_extension([k, l, s], P[i], all_data)
            if s_temp != s:
                curr_list.append([k, l, s])
            if s_temp == 0:
                break
            [k,l,s] = [k_temp, l_temp, s_temp]
        i = i + 1

    # 后向搜索
    prev_list = copy.deepcopy(curr_list) # 深拷贝
    i_temp= len(P)
    i = start_index - 1
    while i >= -1:
        curr_list = []
        s_temp_temp = -1
        for [k, l, s] in prev_list:
            [k_temp, l_temp, s_temp] = backward_extension([k, l, s], P[i] if i != -1 else '$', all_data)
            
            if s_temp == 0 or i == -1:
                if len(curr_list) == 0 and (i+1) < (i_temp+1):
                    i_temp = i
                    match_list.append([k, l, s])
            
            if s_temp != 0 and s_temp != s_temp_temp:
                s_temp_temp = s_temp
                curr_list.append([k, l, s])
        
        if len(curr_list) == 0:
            break

        prev_list = copy.deepcopy(curr_list)

        i = i - 1 

    return match_list
    

    

    

        


