---
layout: mypost
title: linux(debian10)下firefox的flash插件的安装
categories: [linux, flash, firefox]
---

> 前言

&emsp;&emsp;flash虽然被淘汰，但是一些网站仍然使用flash，在linux下的firefox使用flash插件百度上也一堆教程，总体上说就是下载下来插件，解压后放到"~/.mozilla/plugins/"这个目录中，但是好像是自19年开始，flash插件在firefox中被默认关闭，本人在使用manjaro linux的时候，firefox直接不显示它安装了flash插件，在debian10中估计也是如此，不过本人偶然之间从debian wiki找到了解决办法。

> 安装方法

# 命令行
```
apt install browser-plugin-freshplayer-pepperflash
```  
&emsp;&emsp;这条命令在debian jessie(即debian8)以及更新的debian版本中有效，老的debian或者其他发行版，可以使用编译的方法。

# 编译
### 依赖
&emsp;&emsp;先要确保装好了编译所需的依赖：cmake gcc g++ git pkg-config ragel libasound2-dev libssl-dev libglib2.0-dev libconfig-dev libpango1.0-dev libgl1-mesa-dev libevent-dev libgtk2.0-dev libxrandr-dev libxrender-dev libxcursor-dev libv4l-dev libgles2-mesa-dev libavcodec-dev libva-dev libvdpau-dev libdrm-dev libpulse-dev(这些依赖复制自debian wiki，仅供参考)
### 进行编译：
```
cd
git clone https://github.com/i-rinat/freshplayerplugin.git
cd freshplayerplugin
mkdir build
cd build
cmake ..
make
```
### 复制到插件目录
```cp libfreshwrapper-pepperflash.so /usr/lib/mozilla/plugins/```
最后重启一下firefox，就可以测试了。

> 最后一点要说的

* 这个包或者插件实际上是让firefox使用Chromium的flash插件
* 对于其他的发行版或者更老的debian版本，不一定有上述的编译依赖环境，或者说包名不一样，个人认为可以先尝试把这些包安装一下，能装上多少装上多少，全装上自然好，没全装上的话先尝试编译一下，看一下报错信息，再确定应该安装什么包。




