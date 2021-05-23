#include <bits/stdc++.h>
using namespace std;
#define ll long long int
//#include "utility.h"
#define mp make_pair 
const ll mod=1e16+9;
struct display
{
	ll first; ll second; int code;
};

void tokenize(string read,vector<pair<ll,ll>>& kmer,vector<pair<ll,ll>>&hashes,ll k); //form the k-mers and hash them for small k in efficient rolling fashion
void bad_tokenize(string read,vector<pair<ll,ll>>& kmer,vector<pair<ll,ll>>&hashes,ll k); //form the k-mers and hash them for larger k values 
//for larger k values due to overflow issues the text cannot be encoded in rolling fashion

void common1(vector<pair<ll,ll>>& id1,vector<pair<ll,ll>>& id2,
vector<pair<ll,ll>> &common); //find common k-mers between 2 document fingerprints

void get_output(string read,vector<display>& highlight); //display the output with apt highlighting

void seed_and_extend(vector<pair<ll,ll>>& common,vector<pair<ll,ll>>& high1,
vector<pair<ll,ll>>& high2); //used to highlight the overlap using common k-mers

void fill_minimizers(vector<pair<ll,ll>>&hashes,vector<pair<ll,ll>>&min_pos,
ll k,ll w); //extract minimizers from a set of k-meres using monotone stakcs efficiently

string toword(ll lex); //converts a lexicographical code to string

ll murmur64(ll h); //murmur hash algorithm
ll binpow(ll a, ll b,ll m); //binary exponentiation 



set <ll> s1,s2,s3;
bool comp(display d1,display d2)
{
	return (d1.first<d2.first);
	
}
map<ll,ll> fr;

ll func(ll a,ll b)
{
	return a/b;
}
void update(vector<pair<ll,ll>>& hash)
{
	for(auto& x: hash)
	x.first=func(x.first,fr[x.first]);
}
pair<ll,ll> goleft(pair<ll,ll> seed);
pair<ll,ll> goright(pair<ll,ll> seed);
string read1,read2; 
void update();
int main(int argc,char** argv)  //give the 2 input files 
{
    cout<<"Required Density ?"; 
    ll k,w; double den; cin>>den;
    w=2.0/den;
    //w=5;
    fstream in; 
    in.open("text1.txt");
    string line;
    while(getline(in,line))
    read1+=line+"\n";
    if((ll)read1.size()>8)
    k=5;
    else 
    k=2;
    in.close();
    in.open("text2.txt");
    while(getline(in,line))
    read2+=line+"\n";
    in.close();
    
    vector<pair<ll,ll>> kmers1,kmers2,hashes1,hashes2,common,combine,sig1,sig2,high1,high2;
    if(k>=8)
    {
		bad_tokenize(read1,kmers1,hashes1,k);
		bad_tokenize(read2,kmers2,hashes2,k);
	}
	else
	{
		tokenize(read1,kmers1,hashes1,k);
		tokenize(read2,kmers2,hashes2,k);
	}
	
    //common - stores the common representative positions in read 1 and 2
    //high - stores the overlap regions

    fill_minimizers(hashes1,sig1,k,w);
    fill_minimizers(hashes2,sig2,k,w);
    common1(sig1,sig2,common);
    seed_and_extend(common,high1,high2);
    
    vector<display> h1,h2;
    assert(high1.size()==high2.size());
    h1.resize((ll)high1.size()); h2.resize((ll)high1.size());
    
    for(ll i=0;i<(ll)high1.size();i++)
    {
		h1[i].first=high1[i].first; h2[i].first=high2[i].first;
		h1[i].second=high1[i].second; h2[i].second=high2[i].second;
		h1[i].code=i%7; h2[i].code=i%7;
	}
	
	
    sort(h1.begin(),h1.end(),comp);
    sort(h2.begin(),h2.end(),comp);
	
    get_output(read1,h1);
    cout<<"_____________________________________________________________________________"<<endl;
    cout<<endl;
    get_output(read2,h2); cout<<endl;
    cout<<"_____________________________________________________________________________"<<endl;
    cout<<endl;
    
    
    for(auto x: sig1)
    s1.insert(x.first);
    for(auto x:sig2)
    s1.insert(x.first);
    
    
    double jaccard;
    jaccard=s3.size()*1.0/(s1.size());
    cout<<"Jaccard value= "<<jaccard<<endl;

}







ll binpow(ll a, ll b,ll m)
{
    a %= m;
    ll res = 1;
    while (b > 0) {
        if (b & 1)
            res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}
ll murmur64(ll h)
{
	h ^= (h >> 33);
	h *= 0xff51afd7ed558ccdL;
	h ^= (h >> 33);
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= (h >> 33);
	return h;
}
string toword(ll lex)
{
	char x=0;
	if(lex==0)
	return x+"";
	string s;
	while(lex>0)
	{
		int d=lex%128;
		s=(char)d+s;
		lex=lex/128;
	}
	return s;
}
		

void tokenize(string read,vector<pair<ll,ll>>& kmer,vector<pair<ll,ll>>&hashes,ll k)
{
    ll currh=0;
    for(ll i=0;i<k;i++)
    currh=128*currh+read[i];
    ll hpow=binpow(128,k-1,mod);
    kmer.push_back(mp(currh,0));
    for(ll i=k;i<(ll)read.length();i++)
    {
        currh=currh%hpow;
        currh*=128;
        currh+=read[i];
        currh%=mod;
        kmer.push_back(mp(currh,i-k+1));
    }
    for(auto x:kmer)
    hashes.push_back(mp(murmur64(x.first),x.second));
}
void bad_tokenize(string read,vector<pair<ll,ll>>& kmer,vector<pair<ll,ll>>&hashes,ll k)
{
	ll currh;
	for(ll i=0;i+k-1<(ll)read.length();i++)
	{
		currh=0;
		for(ll j=i;j<=i+k-1;j++)
		{
			currh=(128*currh)%mod+read[i];
		}
		kmer.push_back(mp(currh,i));
    }
    for(auto x:kmer)
    {
		hashes.push_back(mp(murmur64(x.first),x.second));
		fr[murmur64(x.first)]++;
	}
}

void common1(vector<pair<ll,ll>>& id1,vector<pair<ll,ll>>& id2,
vector<pair<ll,ll>>& common)
{
    for(auto x:id1)
    {
        for(auto y:id2)
        {
            if(x.first==y.first)
            {
				common.push_back(mp(x.second,y.second));
				s3.insert(x.first);
			}
        }
    }
}
void get_output(string read,vector<display>& highlight) //tested
{
    #define RESET   "\033[0m"
    #define BLUE    "\033[44m"      /* Blue */
    #define mp make_pair
    ll prev,j; prev=0;
	for(ll i=0;i<(ll)highlight.size();i++)
	{
        int col=highlight[i].code; col=41+col;
        string code="\033["+to_string(col)+"m";
        if(prev>highlight[i].first)
        continue;
		for( j=prev;j<highlight[i].first;j++)
		cout<<read[j];
		for( j=max(j,highlight[i].first);j<=highlight[i].second;j++)
        {
		    cout<<code<<read[j];
            cout<<RESET;
        }
		prev=j;
	}

    for(j=prev;j<(ll)read.length();j++)
    cout<<read[j];
}
void seed_and_extend(vector<pair<ll,ll>>& common,
vector<pair<ll,ll>>& high1,vector<pair<ll,ll>>& high2)
{
    for(auto x:common)
    {
        pair<ll,ll> left=goleft(x);
        pair<ll,ll> right=goright(x);
        high1.push_back(mp(left.first,right.first));
        high2.push_back(mp(left.second,right.second));
    }
}
pair<ll,ll> goleft(pair<ll,ll> seed)
{
    ll i1=seed.first; ll i2=seed.second;
    while(read1[i1]==read2[i2])
    {
        i1--; i2--;
        if(i1<0||i2<0)
        break;
    }
    if(i1<0||i2<0||read1[i1]!=read2[i2])
    {i1++; i2++;}
    return mp(i1,i2);
}
pair<ll,ll> goright(pair<ll,ll> seed)
{
    ll i1=seed.first; ll i2=seed.second;
    while(read1[i1]==read2[i2])
    {
        i1++; i2++;
        if(i1>=(ll)read1.length()||i2>=(ll)read2.length())
        break;
    }
    if(i1>=(ll)read1.length()||i2>=(ll)read2.length()||read1[i1]!=read2[i2])
    {i1--; i2--;}
    return mp(i1,i2);
}
void fill_minimizers(vector<pair<ll,ll>>&hashes,vector<pair<ll,ll>>&min_pos,
ll k,ll w) //monotone stacks has been used for better efficiency
{
    ll n=hashes.size();
    ll nse[n],pse[n],i;
    for(i=0;i<n;i++)
    {
        nse[i]=-1;
        pse[i]=-1;
    }
    stack<int> s;
    for(i=0;i<n;i++)
    {
        while((int)s.size()>0&&hashes[s.top()].first>hashes[i].first)
        {
            nse[s.top()]=i;
            s.pop();
        }
        s.push(i);
    }
    s=stack<int>();
    for(i=n-1;i>=0;i--)
    {
        while((int)s.size()>0&&hashes[s.top()].first>hashes[i].first)
        {
            pse[s.top()]=i;
            s.pop();
        }
        s.push(i);
    }
    for(i=0;i<n;i++)
    {
        int prev,next;
        prev=pse[i]; next=nse[i];
        if(prev==-1) prev=0;
        if(next==-1) next=n-1;
        int len=next-prev+1;
        if(len>=w)
        {
			min_pos.push_back(mp(hashes[i].first,i));
		}
        
    }
}


