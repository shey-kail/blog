
### 得到纯文本的方法
因为读入bam奇慢无比，所以最好用读入纯文本的方法，而获得纯文本的方法如下
* processBismarkAln()  
这个是读入bam的函数，与此同时其还有个save.folder的参数，就是指定纯文本生成到何目录
* BismarkCX2methylkit.pl  
这是一个perl脚本，代码如下
```
# BismarkCX2methykit.pl
#NC_000932.1    3   -   0   0   CHH CAT
#NC_000932.1    4   -   0   0   CHH CCA
#NC_000932.1    5   -   0   0   CHH CCC
#NC_000932.1    6   +   0   0   CG  CGA

# chrbase chr base strnd coverage freqC freqT
open (IN,"$ARGV[0]");
open (OUT1,">${ARGV[0]}_CG_methykit.txt");
open (OUT2,">${ARGV[0]}_CHG_methykit.txt");
open (OUT3,">${ARGV[0]}_CHH_methykit.txt");
print OUT1 "chrbase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";
print OUT2 "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";
print OUT3 "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";

while (<IN>) {
   chomp;
@a=split;
$chrbase=$a[0].".".$a[1];
$chr=$a[0];
$base=$a[1];
if ($a[2]=~/-/) { $strand= "R"; 
}
else { $strand ="F";}
$coverage=$a[3]+$a[4];
if ($coverage < 5) {next;}
$freqC= $a[3]/($a[3]+$a[4]) * 100;
$freqT= $a[4]/($a[3]+$a[4]) * 100;
$str=$chrbase."\t".$chr."\t".$base."\t".$strand."\t".$coverage."\t".$freqC."\t".$freqT."\n";
if ($a[5]=~ /CG/) { print OUT1 $str;
} 
elsif($a[5]=~/CHG/) {print OUT2 $str;}
elsif($a[5]=~/CHH/){ print OUT3 $str;}
else{print $_."\n";}
}

```
### reorganize重新提取样本信息
有如下情况，有test1,test2,ctrl1,ctrl2，要比对test1和ctrl2与test2和ctrl1，可以直接用读文件的方法，如下
```
myobj=methRead(file.list=list("test1.txt","ctrl2.txt"),
           sample.id=list("test1","ctrl2"),
           assembly="hg18",
           treatment=c(1,0),
           context="CpG"
           )
meth=unite(myobj, destrand=FALSE)
myDiff=calculateDiffMeth(meth)
all = getMethylDiff(
		myDiff,
		difference = 25,
		qvalue = 0.01,  
		type="all"      
		)
head(all)



myobj=methRead(file.list=list("test2.txt","ctrl1.txt"),
           sample.id=list("test2","ctrl1"),
           assembly="hg18",
           treatment=c(1,0),
           context="CpG"
           )
meth=unite(myobj, destrand=FALSE)
myDiff=calculateDiffMeth(meth)
all = getMethylDiff(
		myDiff,
		difference = 25,
		qvalue = 0.01,  
		type="all"      
		)
head(all)
```
这样读文件读了两遍，挺慢不说还格外占用资源，所以推荐下面的方法
```
### 把所有文件全读入
myobj=methRead(file.list=list("test1.txt","test2.txt","ctrl1.txt","ctrl2.txt"),
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg18",
           treatment=c(1,1,0,0),
           context="CpG",
           mincov = 10
           )

### reorganize重新提取样本信息
myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0))
myobj3=reorganize(myobj,sample.ids=c("test2","ctrl1"),treatment=c(1,0))

meth=unite(myobj2, destrand=FALSE)
myDiff=calculateDiffMeth(meth)
all = getMethylDiff(
		myDiff,
		difference = 25,
		qvalue = 0.01,  
		type="all"      
		)
head(all)
meth=unite(myobj3, destrand=FALSE)
myDiff=calculateDiffMeth(meth)
all = getMethylDiff(
		myDiff,
		difference = 25,
		qvalue = 0.01,  
		type="all"      
		)
head(all)
```
### 多线程
跑整个流程的时候，最慢的环节莫过于读文件和差异甲基化位点的计算（calculateDiffMeth()）,读文件可以靠读纯文本文件来解决，而计算可以通过此函数内的一个多线程相关的参数解决
```
calculateDiffMeth(
			meth,
			mc.cores=2        /*调用两个核心参与运算*/
			)
```
