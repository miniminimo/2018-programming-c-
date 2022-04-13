#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<string>
#include<cmath>
using namespace std;
//------begin--- PART 1-1: Declaration -----------
struct Zone{ //�����U���椧 struct
 int a; //�x�s�Ӻ��椧�a�b�y�СA�Y�ѤW�ӤU�� a �C, a=0,�K,M
 int b; //�x�s�Ӻ��椧��b�y�СA�Y�ѥ��ӥk�� b ��, b=0,�K,N
 int id; //�x�s��s���Aid=0,�K,MN-1
 float r=1.; //�x�s�Ӻ��楻�����ݦs�v�A��l�Ƭ� 1�A
 //�Y�Ӻ��榳���ġA�h�N�� X ��ҳQ�����A��Y�ݯd�F 1-X ��ҡC
 //�խY�t���@���� k����ҧ���ĮĽd��|�i�Υ�����A�h�|�N����
 //k���糧���檺�ݾl�q�]�C�h�@���h���ˮĪG��b�^���ӭ��H
 //�����������ݾl�q�A�H�������C
 bool isM; //�x�s�Y����󥻴�����O�_�Q���ġA�O�h�� 1�A�_�h�� 0
};
struct Virus{ int F; //��ͯf�r�X�{���j��(�Υi�z�Ѭ��ƶq)
 int p; //��ͯf�r�ӱj�פ��o�;��v(%)
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
 int M, //�ϰ��a�b�y�нd��(1,2,�K,M)
 N, //�ϰ��b�y�нd��(1,2,�K,N)
 n_Vir, //��ͯf�r�X�{�������ӼơAĴ�p�ϤG(b)�� n_Vir=3 �N�� 3 ��
        //�f�r�j�סA���O�� 60,40,0 ���v 10%,20%,70%
 T, //�������`����
 s, //���w���f�r�C��������v(%)
 q; //���w���f�r�C���X����v(%)
 Zone *V; //�ʺA�ŧi MN ���פ����� struct �}�C;V[k],k=0,1,�K,MN-1
 Virus *Vir; //�ʺA�ŧi n_Vir ���פ���ͯf�r�}�C;
 float *Q0,//Q0[k],k=0,�K,MN-1 �s���� k �󥻴������`�f�r�ƶq Q0=(Q1+Q2+Q3)*r
 *Q1,//Q1[k]�s k ���W�� Q0 �󥻴������ҭl�ͪ��`�f�r�ƶq Q1=Q0*(1+s)
 *Q2, //Q2[k]�s k ���۾F����W�� Q0 �󥻴��X�����J���f�r�`�q Q2+=�F�� Q0*q
 *Q3;//Q3[k]�x�s���� k �g�L�����ӷs�H���X�{����ͯf�r�ƶq�A���Y�ӱj�� F
 int n_M, //�O��������Ħ��ơA�Y�� 0 �h�N������
 row, //�O����U���Ĥ������a�b�y��
 col, //�O����U���Ĥ������b�y��
 X, //�O����U���Ĥ������v(%)
 R; //�O����U���Ĥ��v�T�d��A�H�ҫ��y�Z����ܡA�Y�� 0 �h�N��ȼv�T���@��
 char filename[50];//�O����J�������ɦW
 int i,j,k;
 //Part1
 readFile(filename,M,N,T,s,q,n_Vir,V,Vir,Q0,Q1,Q2,Q3);
   V=new Zone[M*N];
   for(i=0;i<M*N;i++){
     V[i].id=i;
     V[i].b=i%N;//��b�y�� calculate
     V[i].a=(i-V[i].b)/N;//�a�b�y�� calculate
     V[i].isM=0;//��l������������
   }
   cout<<"Q0[k,0]"<<endl;
   print_colored_2Darray(Q0,M,N);
 //Part2&Part3 �p�G 5 �� n_M=0 �N�� Part2�C�u�n 5 �������@�� n_M!=0 �Y�� Part3
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
            Q0[j]=Q1[j]+Q2[j]+Q3[j];//�[�`�᪺ Q0
        }
        cout<<"Q0[k,"<<k+1<<"]:"<<endl;
        print_colored_2Darray(Q0,M,N);
     }
     else{
        for(i=0;i<n_M;i++){
            cout<<"please input row,col,X,R"<<endl;
            cin>>row>>col>>X>>R;
            residual(X,M,N,V,R,row,col);
            V[row*N+col].isM=1;//�����ĴN�� 1
        }
        growth(Q0,Q1,M,N,s,V);
        expand(Q0,Q2,M,N,s,q,V);
        Original(Q0,Q3,M,N,Vir,n_Vir);
        for(j=0;j<M*N;j++){
            Q0[j]=(Q1[j]+Q2[j]+Q3[j])*V[j].r;//���W�ݦs�v�� Q0
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
 for(k=0;k<M;k++){//amount �@��@��Ū�i�h
    for(i=k*N;i<N*(k+1);i++) {
        file>>Q0[i];
    }
 }
}
void print_colored_2Darray(float *Q0,int &M,int &N){
string col[4]={"\x1b[;32;1m","\x1b[;33;1m","\x1b[;31;1m","\x1b[;35;1m"};
string reset="\x1b[0m";
int i,j,
 interval;//�ƶq���Z
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
//-----begin--- PART 2: Q1,Q2,Q3,Q0(������)-----------
void growth(float *&Q0,float *&Q1,int &M,int &N,int &s,Zone *&V){
 int i,j;
 Q1=new float[M*N];
 for(i=0;i<M*N;i++){
    if(V[i].isM==1)
        Q1[i]=Q0[i];//�p�G�����ĴN�����쥻
    else
        Q1[i]=Q0[i]*(1+s*0.01);//�S���N�� s
 }
}
void expand(float *&Q0,float *&Q2,int &M,int &N,int &s,int &q,Zone *&V){
 int i,j,
 w;//�p��ҫ��y�Z��=1 ���f�r�`�q
 Q2=new float[M*N];
 for(i=0;i<M*N;i++){
    w=0;
    for(j=0;j<M*N;j++){
        if(abs(V[i].a-V[j].a)+abs(V[i].b-V[j].b)<=1&&abs(V[i].a-V[j].a)+abs(V[i].b-V[j].b)>0){//�N�ҫ��y�Z��=1 ���f�r�`�q�[�`�_��(���F�ۤv)
            w+=Q0[j];
        if(V[j].isM==1)
            w=w-Q0[j];//�p�G�d�򤺪����榳�Q���ĴN�
        }
        Q2[i]=w*q*0.01;
    }
 }
}
int probability(Virus *&Vir, int &n_Vir){
 int n,i;
 n=rand()%100+1;
 int *b=new int[n_Vir-1];//�s�ƽu�W�Ʀr�ƶq
 b[0]=Vir[0].p;
 for(i=0;i<n_Vir-1;i++){
    b[i+1]=b[i]+Vir[i+1].p;//�ƽu�W���Ʀr
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
//-----end--- PART 2: Q1,Q2,Q3,Q0(������)-----------
//-----begin--- PART 3:�w���Ī� r-----------
void residual(int X,int M,int N,Zone *&V,int R,int row,int col){
 int i,j,k;
 V[row*N+col].r=1-X*0.01;
 for(j=1;j<=R;j++){
    for(k=0;k<M*N;k++){
        if(abs(row-V[k].a)+abs(col-V[k].b)<=j&&abs(row-V[k].a)+abs(col-V[k].b)>(j-1)){//�p�� r
            V[k].r=V[k].r *(1-(X/pow(2.0,j))*0.01);
        }
    }
 }
}
//-----end--- PART 3:�w���Ī� r-----------
