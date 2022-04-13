#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<string>
#include<cmath>
using namespace std;
//------begin--- PART 1-1: Declaration -----------
struct Zone{ //此為各網格之 struct
 int a; //儲存該網格之縱軸座標，即由上而下第 a 列, a=0,…,M
 int b; //儲存該網格之橫軸座標，即由左而右第 b 行, b=0,…,N
 int id; //儲存其編號，id=0,…,MN-1
 float r=1.; //儲存該網格本期之殘存率，初始化為 1，
 //若該網格有投藥，則將有 X 比例被撲殺，亦即殘留了 1-X 比例。
 //倘若另有一網格 k’其所投放的藥效範圍會波及本網格，則會將網格
 //k’對本網格的殘餘量（每多一單位則殺傷效果減半）拿來乘以
 //本網格原先的殘餘量，以此類推。
 bool isM; //儲存某網格於本期期初是否被投藥，是則為 1，否則為 0
};
struct Virus{ int F; //原生病毒出現之強度(或可理解為數量)
 int p; //原生病毒該強度之發生機率(%)
};
 void readFile(char *filename,int &M,int &N,int &T,int &s,int &q, int &n_Vir,
 Zone *&V,Virus *&Vir,float *&Q0,float *&Q1,float *&Q2,float *&Q3);
 void print_colored_2Darray(float *Q0,int &M,int &N);
 void growth(float *&Q0,float *&Q1,int &M,int &N,int &s,Zone *&V);
 void expand(float *&Q0,float *&Q2,int &M,int &N,int &s,int&q,Zone *&V);
 void Original(float *&Q0,float *&Q3,int &M,int &N,Virus *&Vir,int &n_Vir);
 int probability(Virus *&Vir,int &n_Vir);
 void print_mono_2Darray1(float *Q1,int M,int N);
 void print_mono_2Darray2(float *Q2,int M,int N);
 void print_mono_2Darray3(float *Q3,int M,int N);
 void residual(int X,int M,int N,Zone *&V,int R,int row,int col);
//------end--- PART 1-1: Declaration -----------
int main(){
 int M, //區域縱軸座標範圍(1,2,…,M)
 N, //區域橫軸座標範圍(1,2,…,N)
 n_Vir, //原生病毒出現之種類個數，譬如圖二(b)之 n_Vir=3 代表 3 種
        //病毒強度，分別為 60,40,0 機率 10%,20%,70%
 T, //模擬之總期數
 s, //給定的病毒每期成長比率(%)
 q; //給定的病毒每期擴散比率(%)
 Zone *V; //動態宣告 MN 長度之網格 struct 陣列;V[k],k=0,1,…,MN-1
 Virus *Vir; //動態宣告 n_Vir 長度之原生病毒陣列;
 float *Q0,//Q0[k],k=0,…,MN-1 存網格 k 於本期期末總病毒數量 Q0=(Q1+Q2+Q3)*r
 *Q1,//Q1[k]存 k 之上期 Q0 於本期成長所衍生的總病毒數量 Q1=Q0*(1+s)
 *Q2, //Q2[k]存 k 之相鄰網格上期 Q0 於本期擴散移入之病毒總量 Q2+=鄰格 Q0*q
 *Q3;//Q3[k]儲存網格 k 經過本期而新隨機出現的原生病毒數量，為某個強度 F
 int n_M, //記錄當期投藥次數，若為 0 則代表不投藥
 row, //記錄當下投藥之網格縱軸座標
 col, //記錄當下投藥之網格橫軸座標
 X, //記錄當下投藥之撲殺率(%)
 R; //記錄當下投藥之影響範圍，以曼哈頓距離表示，若為 0 則代表僅影響那一格
 char filename[50];//記錄輸入之測試檔名
 int i,j,k;
 //Part1
 readFile(filename,M,N,T,s,q,n_Vir,V,Vir,Q0,Q1,Q2,Q3);
   V=new Zone[M*N];
   for(i=0;i<M*N;i++){
     V[i].id=i;
     V[i].b=i%N;//橫軸座標 calculate
     V[i].a=(i-V[i].b)/N;//縱軸座標 calculate
     V[i].isM=0;//初始全部都未投藥
   }
   cout<<"Q0[k,0]"<<endl;
   print_colored_2Darray(Q0,M,N);
 //Part2&Part3 如果 5 期 n_M=0 就為 Part2。只要 5 期中有一期 n_M!=0 即為 Part3
  for(k=0;k<T;k++){
    cout<<"please input medicine amounts"<<endl;
    cin>>n_M;
    if(n_M==0){
        cout<<"Q0[k,0]"<<endl;
        print_colored_2Darray(Q0,M,N);
        growth(Q0,Q1,M,N,s,V);
        cout<<"Q1[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray1(Q1,M,N);
        expand(Q0,Q2,M,N,s,q,V);
        cout<<"Q2[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray2(Q2,M,N);
        Original(Q0,Q3,M,N,Vir,n_Vir);
        cout<<"Q3[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray3(Q3,M,N);
        for(j=0;j<M*N;j++){
            Q0[j]=Q1[j]+Q2[j]+Q3[j];//加總後的 Q0
        }
        cout<<"Q0[k,"<<k+1<<"]:"<<endl;
        print_colored_2Darray(Q0,M,N);
     }
     else{
        for(i=0;i<n_M;i++){
            cout<<"please input row,col,X,R"<<endl;
            cin>>row>>col>>X>>R;
            residual(X,M,N,V,R,row,col);
            V[row*N+col].isM=1;//有投藥就為 1
        }
        growth(Q0,Q1,M,N,s,V);
        expand(Q0,Q2,M,N,s,q,V);
        Original(Q0,Q3,M,N,Vir,n_Vir);
        for(j=0;j<M*N;j++){
            Q0[j]=(Q1[j]+Q2[j]+Q3[j])*V[j].r;//乘上殘存率的 Q0
        }
        cout<<"Q1[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray1(Q1,M,N);
        cout<<"Q2[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray2(Q2,M,N);
        cout<<"Q3[k,"<<k+1<<"]:"<<endl;
        print_mono_2Darray3(Q3,M,N);
        cout<<"Q0[k,"<<k+1<<"]:"<<endl;
        print_colored_2Darray(Q0,M,N);
    }
 }
 delete []V;
 delete []Vir;
 delete []Q0;
 delete []Q1;
 delete []Q2;
 delete []Q3;
 return 0;
}
//------begin--- PART 1-2: Read file -----------
void readFile(char *filename,int &M,int &N,int &T,int &s,int &q,int &n_Vir,Zone *&V,Virus *&Vir,float*&Q0,float *&Q1,float *&Q2,float *&Q3){
 int i,j,k;
 fstream file;
 cout << "Please input file name: ";
 cin >>filename;
 file.open(filename, ios::in);
 if(!file)
    cout << "file cannot open";
 file>>M>>N>>T>>s>>q>>n_Vir;
 cout<<"M="<<M<<", N="<<N<<",T= "<<T<<",s="<<s<<"%,q="<<q<<"%,n_Vir="<<n_Vir<<";";
 Vir=new Virus[n_Vir];
 for(i=0;i<n_Vir;i++){
    file>>Vir[i].F;
    file>>Vir[i].p;
 }
 cout<<" F={";
 for(i=0;i<n_Vir;i++){
    cout<<Vir[i].F;
    if(i<n_Vir-1)
        cout<<",";
 }
 cout<<"}, p={";
 for(i=0;i<n_Vir;i++){
    cout<<Vir[i].p;
    if(i<n_Vir-1)
        cout<<",";
 }
 cout<<"};"<<endl;
 Q0=new float[M*N];
 for(k=0;k<M;k++){//amount 一行一行讀進去
    for(i=k*N;i<N*(k+1);i++) {
        file>>Q0[i];
    }
 }
}
void print_colored_2Darray(float *Q0,int &M,int &N){
string col[4]={"\x1b[;32;1m","\x1b[;33;1m","\x1b[;31;1m","\x1b[;35;1m"};
string reset="\x1b[0m";
int i,j,
 interval;//數量間距
char V0='V';
 for(i=0;i<N+1;i++){
    if(i==0)
        printf("%1c",V0);
    else printf("%7.1d",i-1);
 }
 cout<<endl;
 for(i=0;i<M;i++){
    printf("%1.1d",i);
    for(j=i*N;j<N*(i+1);j++){
        if(Q0[j]<33)
            interval=0;
        else if(Q0[j]<66&&Q0[j]>=33)
            interval=1;
        else if(Q0[j]>=66&&Q0[j]<100)
            interval=2;
        else
            interval=3;
        cout <<col[interval] <<setw(7)<<fixed<<setprecision(1)<< Q0[j]<<reset<< flush;
    }
    cout<<endl;
 }
}
//------end--- PART 1-2:Read file -----------
//-----begin--- PART 2: Q1,Q2,Q3,Q0(未投藥)-----------
void growth(float *&Q0,float *&Q1,int &M,int &N,int &s,Zone *&V){
 int i,j;
 Q1=new float[M*N];
 for(i=0;i<M*N;i++){
    if(V[i].isM==1)
        Q1[i]=Q0[i];//如果有投藥就維持原本
    else
        Q1[i]=Q0[i]*(1+s*0.01);//沒有就乘 s
 }
}
void expand(float *&Q0,float *&Q2,int &M,int &N,int &s,int &q,Zone *&V){
 int i,j,
 w;//計算曼哈頓距離=1 的病毒總量
 Q2=new float[M*N];
 for(i=0;i<M*N;i++){
    w=0;
    for(j=0;j<M*N;j++){
        if(abs(V[i].a-V[j].a)+abs(V[i].b-V[j].b)<=1&&abs(V[i].a-V[j].a)+abs(V[i].b-V[j].b)>0){//將曼哈頓距離=1 的病毒總量加總起來(除了自己)
            w+=Q0[j];
        if(V[j].isM==1)
            w=w-Q0[j];//如果範圍內的網格有被投藥就減掉
        }
        Q2[i]=w*q*0.01;
    }
 }
}
int probability(Virus *&Vir, int &n_Vir){
 int n,i;
 n=rand()%100+1;
 int *b=new int[n_Vir-1];//存數線上數字數量
 b[0]=Vir[0].p;
 for(i=0;i<n_Vir-1;i++){
    b[i+1]=b[i]+Vir[i+1].p;//數線上的數字
 }
 if(n<=b[0])
    return Vir[0].F;
 for(i=0;i<n_Vir-1;i++){
    if(n>b[i]&&n<=b[i+1])
        return Vir[i+1].F;
    else
        return Vir[n_Vir-1].F;
 }
}
void Original(float *&Q0,float *&Q3,int &M,int &N,Virus *&Vir ,int &n_Vir){
int i;
 Q3=new float[M*N];
 for(i=0;i<M*N;i++){
    Q3[i]=probability(Vir ,n_Vir);
 }
}
void print_mono_2Darray1(float *Q1,int M,int N){
 int i,j;
 char V0='V';
 for(i=0;i<N+1;i++){
    if(i==0)
        printf("%1c",V0);
    else printf("%7.1d",i-1);
 }
 cout<<endl;
 for(i=0;i<M;i++){
    printf("%1.1d",i);
    for(j=i*N;j<N*(i+1);j++){
        cout<<setw(7)<<fixed<<setprecision(1)<< Q1[j];
    }
    cout<<endl;
 }
}
void print_mono_2Darray2(float *Q2,int M,int N){
 int i,j;
 char V0='V';
 for(i=0;i<N+1;i++){
    if(i==0)
        printf("%1c",V0);
    else
        printf("%7.1d",i-1);
 }
 cout<<endl;
 for(i=0;i<M;i++){
    printf("%1.1d",i);
    for(j=i*N;j<N*(i+1);j++){
        cout<<setw(7)<<fixed<<setprecision(1)<< Q2[j];
    }
    cout<<endl;
 }
}
void print_mono_2Darray3(float *Q3,int M,int N){
 int i,j;
 char V0='V';
 for(i=0;i<N+1;i++){
    if(i==0)
        printf("%1c",V0);
    else
        printf("%7.1d",i-1);
 }
 cout<<endl;
 for(i=0;i<M;i++){
    printf("%1.1d",i);
    for(j=i*N;j<N*(i+1);j++){
        cout<<setw(7)<<fixed<<setprecision(1)<< Q3[j];
    }
    cout<<endl;
 }
}
//-----end--- PART 2: Q1,Q2,Q3,Q0(未投藥)-----------
//-----begin--- PART 3:已投藥的 r-----------
void residual(int X,int M,int N,Zone *&V,int R,int row,int col){
 int i,j,k;
 V[row*N+col].r=1-X*0.01;
 for(j=1;j<=R;j++){
    for(k=0;k<M*N;k++){
        if(abs(row-V[k].a)+abs(col-V[k].b)<=j&&abs(row-V[k].a)+abs(col-V[k].b)>(j-1)){//計算 r
            V[k].r=V[k].r *(1-(X/pow(2.0,j))*0.01);
        }
    }
 }
}
//-----end--- PART 3:已投藥的 r-----------
