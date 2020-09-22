---
layout: mypost
title: debian10安装nvidia闭源驱动和bumblebee纪实
categories: [linux, nvidia]
---

> 前言

&emsp;&emsp;nvidia闭源驱动和bumblebee在linux上的安装是比较麻烦的，关键是会遇见各种各样的坑，今日终于成功，特此纪实，以供参考。

> 第一坑 and 第一步：禁用nouveau

&emsp;&emsp;nouveau是nvidia的开源驱动，很多linux默认自带，它可以让nvidia显卡乖乖的呆在哪里不去闹事，却很难让它高效地去工作（所以预带nouveau是处于稳定性的考虑），但是当你想要安装闭源驱动使你的独显高效工作，nouveau就成了你第一个阻碍，首先你要禁用它：
```
vim /etc/modprobe.d/blacklist.conf
```  
&emsp;&emsp;添加```blacklist nouveau```就可以开机禁用了

> 第二坑 and 第二步：安装nvidia闭源驱动和bumblebee

&emsp;&emsp;事实上，debian的源里面带了nvidia的闭源驱动，用命令行就可以安装，百度上有人说可以使用如下命令：  
```
apt install nvidia-settings
```  
&emsp;&emsp;以前这么安装过，似乎可以，本人因为也同时安装bumblebee，而在debian wiki也说明了，bumblebee和闭源驱动如何一块安装：
```
apt install bumblebee-nvidia primus libgl1-nvidia-glx
```  
&emsp;&emsp;***安装的过程中如果蹦出来什么nvidia-installer-cleanup这个东西，让你yes or no，一定要选no！！！具体yes还是allow什么的我忘了，图也忘了截，以后有机会补上吧***  
&emsp;&emsp;<small><small><small>注意：这是bullseye（也就是debian11）版本之前的版本才可以使用这个，因为在bullseye版本中libgl1-nvidia-glx这个包不存在了，所以使用：```apt install bumblebee-nvidia primus libgl1-nvidia-tesla-glx```来代为安装</small></small></small>

> 第三坑 and 第三步：将用户添加到组

```sudo adduser $USER bumblebee```  
&emsp;&emsp;一条命令解决，这一步其实没啥好说的，如果bumblebee没有的话，自己建一个这样的组就行，然后重启一下。

> 第四坑 and 第四步：安装VitualGL包

&emsp;&emsp;这是最大的一个坑，理论上，你确实安装好了，但是！不妨试试你安装的成果吧！glxgears是nvidia自带的一款测试软件，先```glxgears```测试一下：

![1.png](1.png) 

这是用集成显卡测试的结果，再```optirun glxgears```如果FPS也差不多在这个范围，那么恭喜你你在此坑内，解决办法其实很简单：安装VirtualGL。遗憾的是，debian的源并没有相关的软件包，所以需要现下载，特此提供链接，[点我下载](https://sourceforge.net/projects/virtualgl/files/)，我了amd64的64位的版本和32位的版本，然后  
```dpkg -i software_package.deb```  
就行了

> 第五坑 and 吐槽：CPU占用莫名其妙的很高？？？  

&emsp;&emsp;是的，一直在60%左右，没叫nvidia干活也这样，就是在写本博客的时候发现的，我也很懵逼，累了，以后在解决吧！

> 主要参考文献：  

1. [Debian9 testing 发行版Nvidia驱动+Bumblebee控制显卡切换](https://www.jianshu.com/p/a3161cfa662a)
2. [debian wiki关于bumblebee的描述](https://wiki.debian.org/Bumblebee)
