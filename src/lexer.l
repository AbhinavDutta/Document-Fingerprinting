%option noyywrap
%option yylineno
%{
#include <bits/stdc++.h>
using namespace std;
vector<vector<int>> par,child;
map<string,int> risk;
string lvalue; int curr;
vector<vector<int>> lca;
//constants have risk of -1 by default ie they are safer than any variable
string line;
%}

%x ASSIGN
%x comment

%%
"/*"				BEGIN(comment);
<comment>\n			; 			
<comment>"*"+"/"	BEGIN(INITIAL);
<comment>.			;
[ \t]+				fprintf(yyout," ");
[a-zA-Z]+"="		{   lvalue=yytext; lvalue.pop_back(); curr=-1;
						fprintf(yyout,"R%d=",risk[lvalue]);
						BEGIN(ASSIGN);
					}
<ASSIGN>[a-zA-Z]+	{   if(risk.find(yytext)!=risk.end())
						{
							fprintf(yyout,"R%d",risk[yytext]);
							if(curr==-1)
							curr=risk[yytext];
							else
							curr=lca[curr][risk[yytext]];
						}
					    else
					      fprintf(yyout,"%s",yytext);
				    }
<ASSIGN>;			{
						fprintf(yyout,";");
						risk[lvalue]=curr;
						BEGIN(INITIAL);
					}
<ASSIGN>"(.)*"	fprintf(yyout,"%s","STRING");
<ASSIGN>[0-9]+	fprintf(yyout,"%s","INTEGER");


<ASSIGN>\n			;
<ASSIGN>.			fprintf(yyout,"%s",yytext);

[a-zA-Z]+[a-zA-Z0-9]*		{ 	if(risk.find(yytext)!=risk.end())
					fprintf(yyout,"R%d",risk[yytext]);
					else
					fprintf(yyout,"%s",yytext);
				}
[0-9]+			fprintf(yyout,"%s","INTEGER");
"(.)*"			fprintf(yyout,"%s","STRING");
.				fprintf(yyout,"%s",yytext);
\n				fprintf(yyout,"%s",yytext);

%%
void solver();
int main(int argc,char** argv)
{
	solver();
	extern FILE *yyin, *yyout;
	yyin=fopen(argv[1],"r"); 
	yyout=fopen(argv[2],"w");
	yylex();
}
void dfs(int src,vector<int>& depth,vector<int>&vis)
{
	vis[src]=true;
	for(auto x:child[src])
	{
		if(vis[x])
		continue;
		depth[x]=depth[src]+1;
		dfs(x,depth,vis);
	}
}
void dfs(int src,vector<int>& vis,vector<int>& mark,int sig)
{
	vis[src]=true;
	for(auto x:par[src])
	{
		if(vis[x])
		continue;
		dfs(x,vis,mark,sig);
	}
	mark[src]+=sig;
}


void solver()
{
	fstream read;
	string filename="poset.txt"; 
	read.open(filename);
	int n; read>>n;
	vector<int> depth,vis,mark;
	depth.resize(n+1); vis.resize(n+1); mark.resize(n+1); par.resize(n+1); child.resize(n+1);
	
	//cout<<"Enter number of relations"<<endl;
	int m,i,u,v; read>>m;
	//cout<<"Enter u v if risk u < risk v"<<endl;
	for(i=0;i<m;i++)
	{
		read>>u>>v;
		child[u].push_back(v);
		par[v].push_back(u);
	}
	//cout<<"Enter the least risky element"<<endl;
	int src; read>>src;
	int cnt; string name; int r;
	//cout<<"Enter the number of variable"<<endl;
	read>>cnt;
	for(int i=0;i<cnt;i++)
	{
		read>>name>>r;
		risk[name]=r;
	}

	lca.resize(n+1);
	for(int i=0;i<=n;i++)
	lca[i].resize(n+1);
	
	
	depth[src]=0;
	dfs(src,depth,vis); vis.clear(); vis.resize(n+1);
	
	int lca[n+1][n+1];

	for(int u=0;u<n;u++)
	{
		for(int v=0;v<n;v++)
		{
			mark.clear(); mark.resize(n+1); vis.clear(); vis.resize(n+1);
			dfs(u,vis,mark,1); 
			vis.clear(); vis.resize(n+1);
			dfs(v,vis,mark,2);
			vector<pair<int,int>> node;
			for(int i=0;i<n;i++)
			{
				if(mark[i]==3)
				node.push_back(make_pair(depth[i],i));
			}
			sort(node.begin(),node.end());
			
			if((int)node.size()>0)
			lca[u][v]=(--node.end())->second;
		}
	}
	read.close();

	
}
