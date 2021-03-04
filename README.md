# Shor-Read Aligment Based On The BWT

------

**首先转述原作者的 README 并表达对原作者 [Jwomers](https://github.com/Jwomers) 的感谢**

**源代码链接：https://github.com/Jwomers/burrows_wheeler_alignment**

------

------

#### Burrows-Wheeler Alignment - String Search

------

* This is a very simply quick and dirty implementation of a Burrows-Wheeler Aligner for indexing and sequence alignment.

* It does NOT implement any of the heuristic methods included in the actual BWA algorithm, so implements the basic algorithm which returns 100% accurate results

* It uses the Burrows-Wheeler transform (BWT), Suffix Array (SA), and 2 other auxillary datastructures C and Occ (sometimes called O)

* It uses the D array to prune the tree search and speed up the inexact search algorithm.

* The search is case insensitive.

#### Differences between this code and the real BWA algorithm

------

* This is NOT IN ANY WAY optimised. This has been coded for ease of understanding, ease of reading and with the goal of better understaning the basic BWA algorithm and its datastructures

* This code parses the suffix tree using a recursive depth-first search, while the real BWA uses breadth-first search (and a heap datastructure), and uses this to prioritise the best partial alignments first

* If BWA finds a result with a difference score of z, it only further considers matches with difference z+1, speeding up the process by ignoring worse results and pruning the search space more aggresively

* BWA sets a maximum allowed difference in the first few tens of bases (called the seed) resulting in 2.5x improvement, and very little loss in accuracy. This is not effective on shorter reads (<50bp)

* BWA reduces the required operating memory by storing small fractions of the Occ and SA arrays, and calculating the rest on the fly. This implementation does not do this.

------

------

**首先要说的内容**

使用 python 实现的基于 BWT 的比对算法和实际在 BWA 软件中使用的算法还是很大的区别。

- 如原创作者所述的内容 Differences between this code and the real BWA algorithm。

  具体的细节部分可以参考文献 **Fast and accurate short read alignment with Burrows–Wheeler transform**。

- 使用 python 实现比对算法强化了对算法的理解。实现过程中没有任何**优化**的地方，尽可能清晰表述算法的过程。

  一个循环中仅完成一个步骤/流程。

## 关于代码的说明

整体上分为两个部分

### **1.Class BWA 实现BWT以及精确/非精确匹配**

参考文献1 ***Fast and accurate short read alignment with Burrows–Wheeler transform***
参考文献2 ***Hardware-Acceleration of Short-Read Alignment Based on the Burrows-Wheeler Transform***

`def __init__(self, reference) # 完成BWT相关的数据的计算`

`def OCC(self, char, index, reverse=False) # 获取计算 Occ(a,i) 和 Occ_reverse(a,i) 的值`

  这里值得注意的是：*Occ(a,-1) = 0*

  参数`reverse`判断是  *Occ(a,i)*  还是 *Occ_reverse(a,i)*

`def calculated_D(self, read) # 计算短序列的 D，用于修建遍历树并提升搜寻的速度`

  *D* 的存在可以大幅度缩小便利空间的范围

`def get_D(self, index) # 获取 D 的数值 `

`def InexRecur(self, read, i, z, k, l) #非精确搜索中使用的递归函数`

<br/>
<br/>

按照参考文献1，分别实现的精确匹配和非精确匹配 ，即`def exact_match(self, read)`和`def inexact_match(self, read, z)`。

将两个匹配算法封装在一起，即`def find_match(self, read, num_differences)`。

根据参数`num_differences`判断使用精确还是非精确匹配，显然`num_differences=0`时使用精确匹配。

环境参数`only_InexRecur`的设定可以强制使用非精确匹配算法实现精确匹配，显然此时`z=0`。

`def find_match(self, read, num_differences)`可以认为是**精确匹配 + 非精确匹配**。

<br/>
<br/>

按照参考文献2中提及的算法2，同时实现精确匹配和非精确匹配，即`def all_match(self, read, z)`。

将此算法再次封装，即`def find_match_2(self, read, num_differences)`。

对于该方法来说，可以实现精确匹配就不考虑非精确匹配。

`def find_match_2(self, read, num_differences)`可以认为是**精确/非精确匹配**。

<br/>
<br/>

`def find_match(...)`和``def find_match_2(...)`被封装到一起，即`def match(self, read, num_differences)`。

环境参数`use_find_match`的设定选择使用`def find_match(...)`或者`def find_match_2(...)`。

### **2.Class OCC_COMPRESS 实现相关数据的压缩**

参考文献 ***Hardware-Acceleration of Short-Read Alignment Based on the Burrows-Wheeler Transform***

`def __init__(self, Occ, bwt_str) # 接受 Occ 和 bwt_string 实现数据的压缩`

`def decode(self, char, index) # 对数据进行解码，获取 Occ(a,i) 的值 `

<br/>
<br/>

根据参考文献的表述，数据进行压缩的时候可使用不同的 *entry* 作为 *head*。

环境参数 `use_middle_as_head_number` 设定是否使用每个 *block* 中间的 *entry* 作为 *head*，否的话就是使用第一个 *entry*。

⚠️环境参数`block_number`设定为将 *Occ* 分成的块数，根据文章中的表述为2。

⚠️也可以设定为其他的数值，理论上可以得到更高的压缩率和更小的字符串遍历空间。

⚠️**但是实际上存在问题，随着`block_number`数量的提升，很难做到 *block* 的包含 *entry* 数量的均衡**，目前还没有解决该的问题😢

⚠️**不推荐更改`block_number`的设定值**

## 相关内容的说明

**1.关于 Burrows–Wheeler transform**

参考`./pdf/Burrows–Wheeler transform.pdf`

**2.关于精确匹配算法和非精确匹配算法**

参考`./pdf/Exact&Inexact matching.pdf`

**3.关于 *Occ/O* 辅助数据的压缩**

参考`./pdf/data compression.pdf`

## 疑惑🤔

1.最终得到的匹配结果是所有满足条件的结果的混合，即SNP、indel和精确匹配(此时允许有差异)。

如何将这些匹配结果分开，最主要的是*插入*和*删除*很难以人的角度去辨别。

当差异数相同的时候，需要对所有的匹配结果进行评分和筛选，对目前的我来说是个问题。

2.***Fast and accurate short read alignment with Burrows–Wheeler transform***中认为空串的*SA interval*下限是1，即*k=1*。

但是这样的结论和其他文章有出入，在实际的算法实现中认为 *k=0*。

3.***Hardware-Acceleration of Short-Read Alignment Based on the Burrows-Wheeler Transform***举例*BWT*变换时，对*C*的计算有出入。

算法实现的过程是按照***Fast and accurate short read alignment with Burrows–Wheeler transform***对*C*的定义计算的。

