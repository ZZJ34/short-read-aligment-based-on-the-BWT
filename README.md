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

**1.Class BWA 实现BWT以及精确/非精确匹配**

参考文献 ***Fast and accurate short read alignment with Burrows–Wheeler transform***

**2.Class OCC_COMPRESS 实现相关数据的压缩**

参考文献 ***Hardware-Acceleration of Short-Read Alignment Based on the Burrows-Wheeler Transform***

## 相关内容的说明

**1.关于 Burrows–Wheeler transform**

参考`./pdf/Burrows–Wheeler transform.pdf`

**2.关于精确匹配算法和非精确匹配算法**

参考`./pdf/Exact&Inexact matching.pdf`

**3.关于 *Occ/O* 辅助数据的压缩**

参考`./pdf/data compression.pdf`

## 总结和疑惑🤔