

# Algorithm 1: Sequence retrieval
# Input: Sequence index i≥0; B, O and C defined in the text
# Output: Sequence P and k, the rank of P
def get_seq(index, text_B, text_O, text_C):
    k = index
    p = ''

    while True:
        a = text_B[k]
        k = text_C[a] + text_O[a][k] - 1
        
        if(a == '$'): break

        p = a + p

    return (p, k)

def backward_extension():
    print('hello')

def forward_extension():
    print('hello')


        


