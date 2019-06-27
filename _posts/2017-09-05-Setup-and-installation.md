---
layout: post
title: Setup and installation
---

> ### Windows users: we strongly recommend installing a Unix-like environment such as Cygwin or MinGW to install and use the program. For all intents and purposes, the examples and the commands will be for a Unix environment.

SodaPop is set of command-line tools written in the C++ programming language using the C++11 Standard.

If you are not familiar with a command-line interface, we suggest you read the following section. Otherwise, you may skip directly to the [Installation](#installation).

## Using the command-line

The command-line interface is a powerful way to interact with the computer. In its simplest form, the command-line is a space where you type commands for the computer to execute. On Mac OS X, the command-line is an application called Terminal. It is located in the /Applications/Utilities/ folder. 
The overwhelming majority of command-line programs follow the same syntax. A command can be broken down to three basic components:

  • The utility: also known as the command. In some cases, you can use it without any flag or argument.  
  • The flags: flags are like options. They allow you to modify the behavior of the utility.   
  • The arguments: some utilities take arguments. These are commonly files, but they can be numbers, words or special characters.  

Here is a simple command. We will break it down by component.

```bash
  ls -l Documents/
```  

**ls** is the utility (think of it as a small program). It lists the content of directories.  
**-l** is a flag. It indicates that we want more information than what is provided by default. In fact, think of the l as short for “long”. Flags are most often preceded by a hypen (‘-‘) and consist of single characters. Some flags can also be specified using a word preceded by a double hyphen (‘--‘). Using one or the other will have the same effect. For example, typing

```bash
  man -h
```
or

```bash
  man --help
```

will yield the same response from the terminal. The man utility displays the manual pages for a specific utility. We can call the command with another utility as its argument: 

```bash
  man ls
```  

will display the manual page for **ls**.
The last component of our command is the argument ‘Documents/’. It tells ls that we want to list the content of that directory. And that’s it!  
Besides **ls**, you will need to know how to use a few other basic utilities in order to get started.

```bash
  cd
```  

is the command to change directory. By default, without any flag or argument, **cd** will move up one folder. You can navigate down a folder by giving the name of the folder you wish to move to as an argument. If you are unsure of the folder you are currently in, type

```bash
  pwd
```  

The command stands for *print working directory* and does exactly that.  
If you want to view the contents of a file, say ‘myfile.txt’, type

```bash
  less myfile.txt
```  

This will display it on the screen. You can scroll through the file using the up and down arrows. Pressing Q will quit the utility and bring you back to the command-line.  
If you want to copy a file to a specific location, say ‘Data/’, type

```bash
  cp myfile.txt Data/
```  

Likewise, you can move the file using a similar syntax:

```bash
  mv myfile.txt Data/
```  

Finally, you can remove a file using

```bash
  rm myfile.txt
```  

A helpful feature of the terminal is tab autocompletion. Whenever you type an argument, say a path or a filename, you can type the first few letters and press tab. This will list all the files and folders corresponding to that prefix. If there exists only one, it will autocomplete the argument for you. Getting familiar with this feature will help you to navigate quickly through files and folders in the command-line.

<a name="installation"/>

## Installation

Before you proceed with the installation, you may want to move the compressed repository to a folder of your choosing. To decompress and extract the contents of the downloaded SodaPop repository, open a command-line terminal window, navigate to the directory where the download is located on your computer and execute the following command

```bash
tar –zxvf [name of zip file]
```

You will find executable binaries for all SodaPop tools in the main folder. However, you may want to build the binaries from source. To compile SodaPop on your machine, navigate to the main folder and run

```bash
make
```

followed by

```bash
make clean
```

This will use the makefile to build the binaries **sodapop**, **sodasnap** and **sodasumm**. To install the program to your computer and make it executable globally, run the command

```bash
make install
```

By default, the three components above will be added to /usr/local/bin. You can change this in the makefile by editing the content of the $INSTALLDIR variable. Likewise, any other parameter in the makefile can easily be modified.

If you get the following error when compiling

```bash
  error: unrecognized command line option "-std=c++11"
```  
it is likely that your compiler is out-of-date. You can get a [new version of gcc/g++ here](https://gcc.gnu.org/). Anything from gcc 4.7 onwards will work.

# Great! Now we can [get started]({{ site.baseurl }}{% post_url 2017-09-05-Running-a-basic-simulation %})! 

## or go back to the [home page.]({{ site.baseurl }}{% link index.md %})

