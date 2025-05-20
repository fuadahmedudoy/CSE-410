#include<bits/stdc++.h>
#include <cmath>
#include "matrix.h"
#include "point.h"
#include "Triangle.h"
using namespace std;

point setW(point p){
    point q(p.x,p.y,p.z,p.w);
    q.x/=q.w;
    q.y/=q.w;
    q.z/=q.w;
    q.w/=q.w;
    return q;
}

matrix multiply(matrix A, matrix B){ 
    int row=B.row,col=A.col;  
    matrix temp(row,col);
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
            temp.mat[i][j]=0.0;
            for(int k=0;k<A.col;k++){
                temp.mat[i][j]+=A.mat[i][k]*B.mat[k][j];
            }
        }
    }
    return temp;
}

int main()
{
    ifstream in("scene.txt");
    ofstream o1("stage1.txt");
    ofstream o2("stage2.txt");
    ofstream o3("stage3.txt");`

    double lookx,looky,lookz,eyex,eyey,eyez,upx,upy,upz,fovy,aspectRatio,near,far;    

    stack<matrix> transformations;
    vector<Triangle>triangles;
    matrix M=matrix(4,4),R=matrix(4,4),T=matrix(4,4),P=matrix(4,4);

    in>>eyex>>eyey>>eyez;
    in>>lookx>>looky>>lookz;
    in>>upx>>upy>>upz;
    in>>fovy>>aspectRatio>>near>>far;

    double fovx=fovy *aspectRatio;
    double t=near * tan(fovy *M_PI / 360.0);
    double r=near * tan(fovx *M_PI / 360.0);
    
    double lx,ly,lz,l;
    lx=lookx-eyex;
    ly=looky-eyey;
    lz=lookz-eyez;
    l=sqrt(lx*lx + ly*ly + lz*lz);
    lx/=l,ly/=l,lz/=l;


    double rx,ry,rz,lenr;
    rx=upz*ly - lz*upy;
    ry=upx*lz - lx*upz;
    rz=upy*lx - ly*upx;
    lenr=sqrt(rx*rx + ry*ry + rz*rz);
    rx/=lenr,ry/=lenr,rz/=lenr;

    double ux,uy,uz;
    ux=rz*ly - lz*ry;
    uy=rx*lz - lx*rz;
    uz=ry*lx - ly*rx;

    T.mat[0][3]=-eyex,T.mat[1][3]=-eyey,T.mat[2][3]=-eyez;

    R.mat[0][0]=rx,R.mat[0][1]=ry,R.mat[0][2]=rz;
    R.mat[1][0]=ux,R.mat[1][1]=uy,R.mat[1][2]=uz;
    R.mat[2][0]=-lx,R.mat[2][1]=-ly,R.mat[2][2]=-lz;
    
    P.mat[0][0]=near/r,P.mat[1][1]=near/t,P.mat[2][2]=-(far+near)/(far-near);
    P.mat[2][3]=-(2.0*far*near)/(far-near),P.mat[3][2]=-1.0,P.mat[3][3]=0.0;

    matrix V=multiply(R,T);
    string line;
    int count=0;
    while(1){
        in>>line;
        if(line=="triangle"){
            double x0,y0,z0;
            in>>x0>>y0>>z0;

            double x1,y1,z1;
            in>>x1>>y1>>z1;

            double x2,y2,z2;
            in>>x2>>y2>>z2;

            matrix p0(4,1);
            matrix p1(4,1);
            matrix p2(4,1);

            p0.mat[0][0]=x0,p0.mat[1][0]=y0,p0.mat[2][0]=z0,p0.mat[3][0]=1.0;
            p1.mat[0][0]=x1,p1.mat[1][0]=y1,p1.mat[2][0]=z1,p1.mat[3][0]=1.0;
            p2.mat[0][0]=x2,p2.mat[1][0]=y2,p2.mat[2][0]=z2,p2.mat[3][0]=1.0;

            matrix temp0=multiply(M,p0);
            matrix temp1=multiply(M,p1);
            matrix temp2=multiply(M,p2);

            point transformedp0(temp0.mat[0][0],temp0.mat[1][0],temp0.mat[2][0],temp0.mat[3][0]);
            point transformedp1(temp1.mat[0][0],temp1.mat[1][0],temp1.mat[2][0],temp1.mat[3][0]);
            point transformedp2(temp2.mat[0][0],temp2.mat[1][0],temp2.mat[2][0],temp2.mat[3][0]);

            o1 << transformedp0.x << " " << transformedp0.y << " " << transformedp0.z<<"\n";
            o1 << transformedp1.x << " " << transformedp1.y << " " << transformedp1.z<<"\n";
            o1 << transformedp2.x << " " << transformedp2.y << " " << transformedp2.z<<"\n";
        
            matrix view0=multiply(V,temp0);
            matrix view1=multiply(V,temp1);
            matrix view2=multiply(V,temp2);

            point viewtransformedp0(view0.mat[0][0],view0.mat[1][0],view0.mat[2][0],view0.mat[3][0]);
            point viewtransformedp1(view1.mat[0][0],view1.mat[1][0],view1.mat[2][0],view1.mat[3][0]);
            point viewtransformedp2(view2.mat[0][0],view2.mat[1][0],view2.mat[2][0],view2.mat[3][0]);


            o2 << viewtransformedp0.x << " " << viewtransformedp0.y << " " << viewtransformedp0.z<<"\n";
            o2 << viewtransformedp1.x << " " << viewtransformedp1.y << " " << viewtransformedp1.z<<"\n";
            o2 << viewtransformedp2.x << " " << viewtransformedp2.y << " " << viewtransformedp2.z<<"\n";
            
            matrix project0=multiply(P,view0);
            matrix project1=multiply(P,view1);
            matrix project2=multiply(P,view2);

            point projecttransformedp0(project0.mat[0][0],project0.mat[1][0],project0.mat[2][0],project0.mat[3][0]);
            point projecttransformedp1(project1.mat[0][0],project1.mat[1][0],project1.mat[2][0],project1.mat[3][0]);
            point projecttransformedp2(project2.mat[0][0],project2.mat[1][0],project2.mat[2][0],project2.mat[3][0]);
            
            point pt0=setW(projecttransformedp0);
            point pt1=setW(projecttransformedp1);
            point pt2=setW(projecttransformedp2);

            o3 << pt0.x << " " << pt0.y << " " << pt0.z<<"\n";
            o3 << pt1.x << " " << pt1.y << " " << pt1.z<<"\n";
            o3 << pt2.x << " " << pt2.y << " " << pt2.z<<"\n";        
            
            triangles.emplace_back(pt0,pt1,pt2);

            o1<<"\n";
            o2<<"\n";
            o3<<"\n";

        }
        if(line=="translate")
        {
            double tx,ty,tz;
            in>>tx>>ty>>tz;
            matrix result(4,4);
            result.mat[0][3] = tx;
            result.mat[1][3] = ty;
            result.mat[2][3] = tz;
            M=multiply(M,result);
        
        }
        if(line=="scale")
        {
            double sx,sy,sz;
            in>>sx>>sy>>sz;
            
            matrix result(4, 4);

            result.mat[0][0] = sx;
            result.mat[1][1] = sy;
            result.mat[2][2] = sz;

            M=multiply(M,result);
        }
        if(line=="push")
        {
            transformations.push(M);
        }
        if(line=="rotate")
        {   
            double angle, ax, ay, az;
            in >> angle >> ax >> ay >> az;
            matrix result(4,4);
            double c=cos(angle),s=sin(angle);
            double t=1.0 - c;
            point axis(ax,ay,az,1.0);
            double len = sqrt(ax * ax + ay * ay + az * az);
            axis.x /= len;
            axis.y /= len;
            axis.z /= len;

            point I(1.0, 0.0, 0.0, 1.0);
            point J(0.0, 1.0, 0.0, 1.0);
            point K(0.0, 0.0, 1.0, 1.0);

            point c1 = rodrigues(I, axis, angle);
            point c2 = rodrigues(J, axis, angle);
            point c3 = rodrigues(K, axis, angle);
            
            result.mat[0][0] = c1.x;
            result.mat[0][1] = c2.x;
            result.mat[0][2] = c3.x;
            result.mat[1][0] = c1.y;
            result.mat[1][1] = c2.y;
            result.mat[1][2] = c3.y;
            result.mat[2][0] = c1.z;
            result.mat[2][1] = c2.z;
            result.mat[2][2] = c3.z;
            M=multiply(M,result);
        }
        if(line=="pop"){
            M=transformations.top();
            transformations.pop();
        }
        if(line=="end"){
            break;
        }
    }
    in.close();
    o1.close();
    o2.close();
    o3.close();
}