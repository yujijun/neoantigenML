# This is a document for work recode.
0. 如何手动创建一个R包？
参考印象笔记：“如何动手创建一个R包”

1. 如何去加载自己编写的函数，在不同的function 之间？
A： 用source + 相对路径进行加载

2. 如何加载别人的包？
A： 直接library就行

3. 关于如何创建R包并提交到github，可以什么文件？
A: https://www.jianshu.com/p/363e867de465

4. 如何进行序列编码：
DNA 序列
4.1 顺序编码
4.2 one-hot 编码
4.3 将DNA序列视为语言，进行K-mer计算，然后按照词汇分析的方法进行分析
https://blog.csdn.net/weixin_44022515/article/details/103867115
https://www.cxyzjd.com/article/qq_36333576/106497423

蛋白质序列
https://zhuanlan.zhihu.com/p/34562954 蛋白质序列进行独热编码

4.1 R 中如何执行one hot 编码：
https://datatricks.co.uk/one-hot-encoding-in-r-three-simple-methods


这部分可以等到模型训练好之后，接着进行motif 的开发：
4. 那些包可以用在peptide 聚类 和 motif 发现？
A: 
(1) BiocManager::install("memes") 
(1.1) seqinr： Biological Sequences Retrieval and Analysis

(2) 计算一些peptide的生物学属性和指标 https://cran.r-project.org/web/packages/Peptides/Peptides.pdf

(3) Biostrings: Efficient manipulation of biological strings; 可以用来序列读取，计算PWM

(4) 计算motif可以参考：https://bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html
https://www.bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack_HTML.html#Examples_of_using_motifStack
(4.1) 展示motif 可以参考：seqLogo：https://bioconductor.org/packages/release/bioc/vignettes/seqLogo/inst/doc/seqLogo.html

(5) 计算距离用Rdist

rdist 常用来计算曼哈顿距离

思考几种计算距离的方式：
这些其实没有涉及到距离的计算：
1. PFM(Position Frequency Matrix):/PPM(position probability matrix)
2. PWM(position weight matrix)； #从motif进入
3. 通过PWM对序列进行最终计算得分，看是否是一个确定的motif。确定最大的motif


The TCRdist distance between two TCRs is defined to be the similarity-weighted mismatch distance, The mismatch distance is defined based on the BLOSUM6237 substitution matrix. - 计算出一个距离矩阵。


(6) 计算聚类用CLIPH
a collection of all continuous 2-mers, 3-mers, 4-mers and 5-mers + fisher exact test

关于聚类计算的可信度：引入噪音 + 计算聚类纯度。如果超过90% 则说明计算比较好。

6. 关于R function 上面的参数？
#' @title Poisson vector.
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer.
#' @param x A non-negative integer. or vector.
#' @return A numeric vector of log density.
#' @export
#' @import stats
#' @importFrom stats dpois
#' @examples
#' poiVec(5)
#' poiVec(6)
poiVec <- function(x) {
  dpois(x , x)
}

1）注释的每一行开头用#'或者##'开始。注释基本上都是完整的句子，请用实心句号结束。

2）注释的内容和格式，每一个项目和其内容之间用空格隔开：

@title 定义函数的标题，这个标题里的词可以搜索到

@description 是函数的描述，如果描述内容很长，则可以写在@details里。

@param 声明函数的参数，每个参数都需要声明，如果用一行同时声明2个输入类型一致的
参数，则用@param x,y 来表示，参数之间只用逗号隔开，不能有空格。

@return 非必须，声明返回的内容和类型，如果没有返回就不需要加这个项，否则会报错。

@export 非必须，如果需要输出该函数到NAMESPACE里，则必须加上这个项

@import 和@importFrom 虽然不是必须的，但是非常需要注意，如果你的函数里用到了其他包的函数，最好声明一下，当然基本的函数可以不用声明。这些内容会写到NAMESPACE文件里。如果不声明的话，在check包的时候，就会因为运行不了函数的例子而报错。

而且不建议在函数里加library()或者require()这种命令，也不建议在函数里使用包名::函数名这种命令，最好通通写在import和importFrom里，这样在生成rd文件的时候就被roxygen2自动输出到NAMESPACE文件中，方便维护和修改。

@examples 也不是必须的，但是一旦你写了，就要保证里面的脚本能够运行。所以一定要声明好import的包。还有如果example里用到了其他的包，需要写require(包名)，否则check的时候会找不到这个包。

8. How to add R script header template?
https://timfarewell.co.uk/my-r-script-header-template/
https://www.rdataguy.com/2019/11/lesson-7-creating-rscript-header_11.html


