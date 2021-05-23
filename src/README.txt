
USAGE INSTRUCTIONS

TEXTUAL INPUT
Step 1: Enter the 2 files in text1.txt and text2.txt
Step 2: Run the main.cpp ($g++ main.cpp $./a.out)

CODE INPUT
Step 1: Enter the 2 files in in1.txt and in2.txt
Step 2: Enter the poset - related information in poset.txt
Step 3: Mask the files and write the masked files in text1.txt and text2.txt
		$flex lexer.l
		$g++ lex.yy.c
		$./a.out in1.txt text1.txt
		The above code takes the code from in1.txt and writes the masked code to text1.txt
		Similar process should be followed to mask in2.txt and write to text2.txt
Step 4: (Optional - If you want to detect memory issues in the masked code in text1.txt )
		$flex memory.l
		$g++ lex.yy.c
		$./a.out text1.txt
Step 5: Run main.cpp to highlight the matching text and calculate the similarity metrics like jaccard index.
		$g++ main.cpp
		$./a.out
		
____________________________________________________________________________________________

The purpose of the files are listed below :

1)in1.txt,in1.txt - If the input is a program code then enter the input here
Finally either the codes of in1.txt and in2.txt will be compared for similarity or detection
of memory issues.

2)poset.txt - Enter the information about the poset in this file. The format is follows:
number of categories
number of egdes in the hasse diagram of the poset
edges in (u,v) if (u<=v in the poset 
smallest element in the poset
number of variables in the program
the variables and their associated category 
Example (A simple poset with 2 categories where category0<=category1)
2
1
0 1 
0
3
alpha 0
beta 0
gamma 1

3)memory.l - This is a lexer that detects memory issues like memory leak and dangling pointer.
This assumes that all the lines in the program will be executed and thus might be inaccurate in practice

4)lexer.l - This is a lexer that does two things :
			a) Remove whitespace and comments 
			b)Implement the masking scheme and write the masked code in output file
5)main.cpp - This contains the implementation of the minimizers and the frequency weighted minimizers.
			It does the following :
			1) Computes the k-mers and hash them (murmur hash algorithm is used)
			2) Choses a subset of the k-mers as representatives of the document. This smaller set of k-mers
			is called the document fingerprint. There are many schemes for selecting a subset of k-mers. Two of
			them have been explored - minimizers (used in MOSS) , frequency weighted minimizers (recently discovered in 2018
			ans used in a genetic analysis tool called Winnowmap)
			3) Common k-mers from the 2 files are recorded
			4) Seed - and extend to highlight overlaps.

6)text1.txt,text2.txt - If the input is a textual document then enter it here. Actually, program codes after masking are also 
written these files.
____________________________________________________________________________________________________________________________
FEATURES

1) There are 2 parameters - k and w. k denotes the k-mer length and dictates the selectivity-specificity tradeoff.
Too large or too small values of k will give no match or many random matches. By default I have set k=8 . However for
larger files (1000+ characters) k=15 is recommended .Also while dealing with very large data-sets k=31 seems spt.
w on the other hand primarily affects how much computational resources will be used. This has a cost-accuracy tradeoff. While dealing with
large set of large file selecting w=100-1000 seems apt because in such cases we only want to detect large matches. 

2) Frequency weighted minimizers- If the text contains lot of large repeated words then it is recommended to use frequency weighted
minimizers. Though there is objective way of measuring how much better this scheme is over normal minimizers due to lack of datasets
however it is ecxpected to do better in such cases. By default, the program uses normal minimizers but there is an in-built update function
to re-vaule the minimizers based on frequency and that can be used in this case

_____________________________________________________________________________________________________________________________-
