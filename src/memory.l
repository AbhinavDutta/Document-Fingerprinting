%option noyywrap
%option yylineno
%{
#include <bits/stdc++.h>
using namespace std;
map<string,int> address; //map variables to addr.
map<string,int> category; //pointers have category 1 and others 2
set<string> location[1000]; //all the variables assc. with location[i]
map<string,int> dang_mark,leak_mark;
string lvalue; string line;
int cnt=0;
%}

%x ASSIGN
%x DEALLOC

%%			

"free"|"delete"		{BEGIN(DEALLOC);}
<DEALLOC>"("		;
<DEALLOC>[a-zA-Z0-9]+	{
						int add=address[yytext];
						if(add>=0){
						fprintf(yyout,"Potential Dangling Pointers at line %d : ",yylineno);
						for(auto x:location[add])
						{
							fprintf(yyout,"%s ",x.c_str());
						}
						fprintf(yyout,"\n");}
						address[yytext]=-1;
					}
<DEALLOC>")"		BEGIN(INITIAL);
[a-zA-Z0-9]+"=malloc"	{
						lvalue=yytext; 
						for(int i=0;i<7;i++)
						lvalue.pop_back();
						if(address[lvalue]>=0)
						location[address[lvalue]].erase(lvalue);
						cnt++; address[lvalue]=cnt;
						location[cnt].insert(lvalue);
					}
[a-zA-Z0-9]+"=NULL"	{ 	lvalue=yytext; 
						for(int i=0;i<5;i++)
						lvalue.pop_back();
						if(address[lvalue]>=0)
						location[address[lvalue]].erase(lvalue);
						address[lvalue]=-1;
					}
						
[a-zA-Z0-9]+"="		{   lvalue=yytext; lvalue.pop_back();
						if(category[lvalue]==1)
						{
						int add=address[lvalue];
						if(add!=-1&&(int)location[add].size()<=1)
						{
							fprintf(yyout,"Potential Memory Leak at line %d due to %s\n",yylineno,lvalue.c_str());
						}
						}
						//fprintf(yyout,"%s\n",lvalue.c_str());
						BEGIN(ASSIGN);
					}
<ASSIGN>[a-zA-Z0-9]+	{	int add=address[yytext];
						if(address[lvalue]>=0)
						location[address[lvalue]].erase(lvalue);
						address[lvalue]=add;
						location[add].insert(lvalue);
					}
		
<ASSIGN>;			{
						BEGIN(INITIAL);
					}
<ASSIGN>["](.)*["]		{ 	
						int add=address[lvalue];
						if(add!=-1&&(int)location[add].size()<=1)
						{
							fprintf(yyout,"Potential Memory Leak at line %d due to %s\n",yylineno,lvalue.c_str());
						}
						cnt++;
						if(address[lvalue]>=0)
						location[address[lvalue]].erase(lvalue);
						address[lvalue]=cnt; 
						location[cnt].insert(lvalue);
					}
<ASSIGN>\n		
<ASSIGN>.			;

.					;
\n					;

%%
void solver();
int main(int argc,char** argv)
{
	solver();
	extern FILE *yyin, *yyout;
	yyin=fopen(argv[1],"r"); 
	yyout=fopen(argv[2],"w");
	yylex();
	cout<<"success"<<endl;
}
