### 安装
```
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("methyLkit")
```
### 使用
methyLkit支持两种文件的输入
* 纯文本格式
* bam


纯文本格式的文件，其内容遵循如下格式  
```
chrBase chr base strand coverage freqC freqT 
chr21.9764539 chr21 9764539 R 12 25.00 75.00 
chr21.9764513 chr21 9764513 R 12 0.00 100.00 
chr21.9820622 chr21 9820622 F 13 0.00 100.00 
chr21.9837545 chr21 9837545 F 11 0.00 100.00 
chr21.9849022 chr21 9849022 F 124 72.58 27.42 
chr21.9853326 chr21 9853326 F 17 70.59 29.41 
```
每一行是一个甲基化位点，coverage 代表覆盖这个位点的reads数，freqC 代表甲基化C的比例，freqT 代表非甲基化C的比例。这种纯文本格式内容非常直观，文件大小相比bam 文件小很多，读取的速度更快。
#### 读取
* 纯文本的文件的读取
```
library(methylkit)

//file.list这里存放的是每个样品的纯文本文件的路径
file.list = list("sample_A1.txt",
		"sample_A2.txt",
		"sample_B1.txt",
		"sample_B2.txt"
		)

//id是给每个文本文件赋予一个对应的id
id = list(
		"sample_A1",
		"sample_A2",
		"sample_B1",
		"sample_B2"
		)


myobj = methRead(
		location = file.list,
		sample_id = id,
		assembly = "hg19",         //这个参数是描述基因组的一个字符串，可以是任何值
		treatment = c(1,1,0,0),    //在这里，1代表实验组，0代表对照组
		context="CpG"
		)
```
* bam文件  
其实bam文件的读取跟上面区别不大，就是用了个别的函数
```
library(methylkit)

//file.list的作用同上
file.list = list("sample_A1.bam",
		"sample_A2.bam",
		"sample_B1.bam",
		"sample_B2.bam"
		)

//也同上
id = list(
		"sample_A1",
		"sample_A2",
		"sample_B1",
		"sample_B2"
		)


myobj = processBismarkAln(
		location = file.list,
		sample_id = id,
		assembly = "hg19",         //这个参数是描述基因组的一个字符串，可以是任何值
		treatment = c(1,1,0,0),    //在这里，1代表实验组，0代表对照组
		context="CpG",
		save.folder = getwd()      //暂时不清楚为啥加这东西
		)
```
#### 合并所有情况，得到所有样本的甲基化表达谱
很简单，只用一个函数  
`meth=unite(myobj, destrand=FALSE)`  
所得的meth可以这样看  
`head(meth)`  
实际上就是上面样品文件所含信息的合并
#### 差异化位点分析
`myDiff=calculateDiffMeth(meth)`  
根据甲基化C是变多了还是变少了，可以将差异甲基化的结果分为两大类：
* **hypermethylated**
* **hypomethylated** 

**hypermethylated**表示相比对照组，实验组中的甲基化C更多；**hypomethylated**则相反，表示实验组中的甲基化C比对照组中少。  
采用getMethylDiff函数提取差异分析的结果，用法如下
```
all = getMethylDiff(
		myDiff,
		difference = 25,  //这个东西是差异阈值，只有差异大于阈值时，才会保留
		qvalue = 0.01,    //q值，小于它时才保留
		type="all"        //这个参数的取值范围为all,hypomethylated,hypermethylated，用于告诉R你想要显示那部分的结果
		)
```
#### 差异甲基化区域的分析
在methylKit中，它的差异分析总是针对合并后的甲基化表达谱，如果你的甲基化表达谱每一行是一个甲基化位点，那么差异分析的结果就是差异甲基化位点；如果你的表达谱每一行是一个甲基化区域，那么差异分析的结果就是差异甲基化区域。上面的例子都是针对差异甲基化位点的，下面看下差异甲基化区域的分析。

首先遇到的问题就是甲基化区域如何界定，在methylKit中，按照滑动窗口的方式定义甲基化区域，默认窗口大小为10000 bp ，步长为1000bp,通过**tileMethylCounts**函数实现。

完整的差异甲基化区域分析的代码如下：
```
//以下步骤是读取完样品文件以后做的事
regions = tileMethylCounts(myobj,win.size=1000,step.size=1000)
meth = unite(regions,destrand=FALSE)
myDiff = calculateDiffMeth(meth)
all = getMethylDiff(
		myDiff,
		difference = 25,  
		qvalue = 0.01,    
		type="all"        
		)
head(all)
```
