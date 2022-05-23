# Studying with debug
I can't debug successfully. But , I also don't know what went wrong.   

When I have enough experience, I will try to revise.

# Failed to debug

# Succeeded to debug

## ./CUDA/matrix_multiplication.cu 
**What did I learn?**    
&emsp;&emsp;I should need a `tool` to develop a flowchart before I write a `coding scheme`.  
&emsp;&emsp;I don't consider logic in code as a diffculty ,and I rarely fail in logic when I write a code in reality. But, I usually failed in a tiny thing what is spelling mistake or using other variable . I have correct logic in my mind , but my hand or my `muscle memory` give the code a fatal error. It's afflictive !  
&emsp;&emsp;So I think that I need a coding scheme to help me to check code to resisting the bad effects of muscle memory.   
&emsp;&emsp;I will use `Visio` to do this. And I hope I can find better tool.

## ./CT/Split_Bregman 2022/5/19
&emsp;&emsp;The program took me a week. On roads to coding-success, it's afflictive 😟 On roads to success in the experiment which gived me good results, it breaks my faith 😳. Through this exprience, I have learned a lot. Then, I will record my thoughts.  
**What did I learn?**  
* ___`Functional Programming`/`Structuerd Programming(SP)` is the better programming paradigm than `Procedural Programming`.___<br>
&emsp;&emsp;Some time before that, I don't like Functional Programming. Because I think Procedural Pragramming makes me feel comfortable when I write a code, and Functional Programming needs too many parameters to write into Function particularly when you need nested function calls.<br>
&emsp;&emsp;However, It's a miserable thing that review a code using Procedural Pragramming and having over 509 lines. When I debug a code using Functional Pragramming, the review is so easy.**By using Functional Pragramming, the work of review is clear: check the Spelling mistake and logic error in this function; check the parameters about both type and definition; check the value of return results. Through focusing on these points, you can review your code clearly and fastly.** I was impressed (smiling_face_with_tear).  
* ___“Easy to read, review, and extend” is the key to write a good code. If you not need to consider high-performance , many techniques aren't important!___

## ./FILE_stream/imgData.dat 2022/5/23
&emsp;&emsp;My teacher(Yining Zhu) gived us a picture of the Raw Data type a week ago, which I've never used data of this type. 
&emsp;&emsp;The characteristic of  this picture: 
- Gray Image;
- [width, height] = [4096, 2048];
- Single-precision floating point type.  

&emsp;&emsp;He also tell us that we can use `ImageJ` to read this picture.But in the beginning, I can't sucess to read this picture, through I think I have know all information about this problem.  
&emsp;&emsp;I also have tried to use `<cstdio>` to solve this problem. But, all the results of gray value of each pixel  is zero. I'm very confused!  
**What did I learn?**  
* ___`Real` also is a Single-precision floating point type.___
* ___HOW to use `<cstdio>` to read file.___  
    - The basic grammer of `<cstdio>`.  
    - Along with using C/C++, the more I feel that C++ has many advantages over C in grammer. For example, when I an using  `size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)`, I usually consider that `size_t size` or `size_t nmemb` is redundant. In reality, we usually make one of these equal to 1. 😅  
* ___The range of the gray value of Single-precision floating point type image is [0,1].___  
 &emsp;&emsp;So when you want to get gray image, you should multiply the gray value by `255`.  
 
