#This is an neoantigens filtering repository by machine learning method

** Description about function:**
1. Main function for whole project was in main.r function;
2. All preprocess about datasets were inclueded into dataset.r;
3. All functions were put into *.R;
4. All supplementary functions used in *.R were included into base.r;
5. All Needed packaged was included in InstallPackage.R, User could install all packages by running this script.
** How to change code? (will be deleted when we would like to publish)**

1. 打开dataset.r 整理数据，将需要加载到环境中的数据，用use_data 函数进行加载；
2. 打开*.R 编写运行函数；
3. 将数据 和 *.R 函数 在mian.r 函数中运行【不要将load_all 函数写入main.r 函数】；
4. document() 将更新的函数进行编译；
5. 将更改提交到github.

Attention:
1. 注意每次更改完函数之后一定要用document()提交编译.
2. 将main.r 函数不要放到 ./R folder 下面，编译过程会直接提交运行,将其单独放到main_running folder.
3. 如果涉及到自己的写的函数之间的调用，注意使用source函数，写清楚相对路径.


关于结果文件的统计：用main_running/ folder 底下的deredundancy.R script，去处理IEDB 中Tcell_full_v3 文件，主要包含一下几个步骤：
 ---------------
1. 列名规整，将文件的列名，用janitor 函数进行规整；
2. 文件删除descirption\qualitative_measure\allele_name 中任何一列包含NA数据的term；删除非线性肽段 和 肽段中包含[+]新的term；
3. 肽段长度选择： 计算肽段长度，并且筛选出肽段长度为9的肽段，然后留下最重要的几列信息，description,allele_name,qualitative_measure,organism_name,name,pep_length
4. 肽段去重复：对于完全相同的肽段，如果有一个证明为阳性，则将阳性结果留下，如果有超过一个阳性，判断MHC 亚型，将所有不同亚型的阳性结果留下；如果没有阳性结果，将所有不同MHC亚型的阴性结果留下。【注意：此地没有将同样长度，包含同样结合位点的肽段冗余情况考虑进去，如：AAANIIRTL 和 AANIIRTN 等】
5. 添加判断的judge标签：将有一个阳性的peptide标注为1，没有任何阳性标签的标注为0 
6. 肽段distinct：将相同肽段不同MHC情况的term 只留下其中一个term【如：AAAAQQIQV 在HLA-A* 02:01 和 HLA-B* 07:02中都为 Negative 情况，我们只留一个进行进一步处理】
7. 最后筛选出21537 个 unique peptides for further analysis.  结果保存在：/Users/yujijun/Documents/work/4_AntigenML/NeoantigenML/result/yjj/Tcell_deredundancy_all.distinct.RData 文件中; hosthuman 和 hostMus 主要保存在Tcell_deredundancy_hosthuman.RData 和 Tcell_deredundancy_hostMus.RData 文件中
8. 通过blast的方式进行肽段去冗余，具体过程见下图。
------------

具体示意图如下：
![](https://github.com/yujijun/neoantigenML/blob/main/image/deredundancy.jpg)

Reference:
1. please refer to https://mlr3book.mlr-org.com/ for detail usage about mlr3. 
