#This is an neoantigens filtering repository by machine learning method

** Description about function:**
1. Main function for whole project was in main.r function;
2. All preprocess about datasets were inclueded into dataset.r;
3. All functions were put into *.R;
4. All supplementary functions used in *.R were included into base.r;

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

Reference:
1. please refer to https://mlr3book.mlr-org.com/ for detail usage about mlr3. 
