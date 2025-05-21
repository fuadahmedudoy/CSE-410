#include<bits/stdc++.h>
#include <cmath>
#include "matrix.h"
#include "point.h"
#include "Triangle.h"
#include "bitmap_image.hpp"
using namespace std;



point setW(point p){
    point q(p.x,p.y,p.z,p.w);
    q.x/=q.w;
    q.y/=q.w;
    q.z/=q.w;
    q.w/=q.w;
    return q;
}
bool EQ(double x, double y){
    return fabs((x) - (y) ) < __DBL_EPSILON__; 
}

matrix multiply(matrix A, matrix B){ 
    int row=A.row,col=B.col;  
    matrix temp(row,col);
    for(int i=0;i<A.row;i++){
        for(int j=0;j<B.col;j++){
            temp.mat[i][j]=0.0;
            for(int k=0;k<A.col;k++){
                temp.mat[i][j]+=A.mat[i][k]*B.mat[k][j];
            }
        }
    }
    return temp;
}
double getYMinPoint(Triangle t)
{
    double min1=min(t.points[0].y,t.points[1].y);
    double ans=min(min1,t.points[2].y);
    return ans;
}
double getXMinPoint(Triangle t)
{
    double min1=min(t.points[0].x,t.points[1].x);
    double ans=min(min1,t.points[2].x);
    return ans;
}
double getZMinPoint(Triangle t)
{
    double min1=min(t.points[0].z,t.points[1].z);
    double ans=min(min1,t.points[2].z);
    return ans;
}
double getYMaxPoint(Triangle t)
{
    double max1=max(t.points[0].y,t.points[1].y);
    double ans=max(max1,t.points[2].y);
    return ans;
}
double getXMaxPoint(Triangle t)
{
    double max1=max(t.points[0].x,t.points[1].x);
    double ans=max(max1,t.points[2].x);
    return ans;
}
double getZMaxPoint(Triangle t)
{
    double max1=max(t.points[0].z,t.points[1].z);
    double ans=max(max1,t.points[2].z);
    return ans;
}

int main()
{
    ifstream in("scene.txt");
    ofstream o1("stage1.txt");
    ofstream o2("stage2.txt");
    ofstream o3("stage3.txt");

    double lookx,looky,lookz,eyex,eyey,eyez,upx,upy,upz,fovy,aspectRatio,near,far;    

    stack<matrix> transformations;
    vector<Triangle>triangles;
    matrix M=matrix(4,4),R=matrix(4,4),T=matrix(4,4),P=matrix(4,4);

    in>>eyex>>eyey>>eyez;
    in>>lookx>>looky>>lookz;
    in>>upx>>upy>>upz;
    in>>fovy>>aspectRatio>>near>>far;

    double fovx=fovy * aspectRatio;
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
    ux=ry*lz - lz*ry;
    uy=rz*lx - lz*rx;
    uz=rx*ly - lx*ry;

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

            o1<<fixed<<setprecision(7) << transformedp0.x << " " << transformedp0.y << " " << transformedp0.z<<"\n";
            o1<<fixed<<setprecision(7) << transformedp1.x << " " << transformedp1.y << " " << transformedp1.z<<"\n";
            o1 <<fixed<<setprecision(7)<< transformedp2.x << " " << transformedp2.y << " " << transformedp2.z<<"\n";
        
            matrix view0=multiply(V,temp0);
            matrix view1=multiply(V,temp1);
            matrix view2=multiply(V,temp2);

            point viewtransformedp0(view0.mat[0][0],view0.mat[1][0],view0.mat[2][0],view0.mat[3][0]);
            point viewtransformedp1(view1.mat[0][0],view1.mat[1][0],view1.mat[2][0],view1.mat[3][0]);
            point viewtransformedp2(view2.mat[0][0],view2.mat[1][0],view2.mat[2][0],view2.mat[3][0]);


            o2 <<fixed<<setprecision(7)<< viewtransformedp0.x << " " << viewtransformedp0.y << " " << viewtransformedp0.z<<"\n";
            o2 <<fixed<<setprecision(7)<< viewtransformedp1.x << " " << viewtransformedp1.y << " " << viewtransformedp1.z<<"\n";
            o2 <<fixed<<setprecision(7)<< viewtransformedp2.x << " " << viewtransformedp2.y << " " << viewtransformedp2.z<<"\n";
            
            matrix project0=multiply(P,view0);
            matrix project1=multiply(P,view1);
            matrix project2=multiply(P,view2);

            point projecttransformedp0(project0.mat[0][0],project0.mat[1][0],project0.mat[2][0],project0.mat[3][0]);
            point projecttransformedp1(project1.mat[0][0],project1.mat[1][0],project1.mat[2][0],project1.mat[3][0]);
            point projecttransformedp2(project2.mat[0][0],project2.mat[1][0],project2.mat[2][0],project2.mat[3][0]);
            
            point pt0=setW(projecttransformedp0);
            point pt1=setW(projecttransformedp1);
            point pt2=setW(projecttransformedp2);

            o3 <<fixed<<setprecision(7)<< pt0.x << " " << pt0.y << " " << pt0.z<<"\n";
            o3 <<fixed<<setprecision(7)<< pt1.x << " " << pt1.y << " " << pt1.z<<"\n";
            o3 <<fixed<<setprecision(7)<< pt2.x << " " << pt2.y << " " << pt2.z<<"\n";        
            
            triangles.emplace_back(pt0,pt1,pt2);
            // cout<<triangles.size()<<endl;

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
    in.open("config.txt");
    ofstream o4("z_buffer.txt");
    int width ,height;
    double leftLimit,rightLimit,topLimit,bottomLimit,zMin,zMax;
    in>>width>>height;
    in>>leftLimit;
    in>>bottomLimit;
    in>>zMin>>zMax;
    rightLimit=-1.0*leftLimit;
    topLimit=-1.0*bottomLimit;

    double dx = (rightLimit - leftLimit) / width,dy = (topLimit - bottomLimit) / height;
    double leftX = leftLimit + dx / 2.0,rightX = rightLimit - dx / 2.0;
    double topY = topLimit - dy / 2.0,bottomY = bottomLimit + dy / 2.0;

    double **zbuffer = new double *[height];
    for(int i=0;i<height;i++){
        zbuffer[i]=new double[width];
        for(int j=0;j<width;j++){
            zbuffer[i][j]=zMax;
        }
    }
    //import first
    bitmap_image image(width, height);
    image.set_all_channels(0, 0, 0);
    //cout<<"befr loop\n";
    //cout<<triangles.size();    
    for(Triangle t: triangles)
    {
        //cout<<t.points[0].x<<" "<<t.points[1].y<<" "<<t.points[2].z<<endl;;
        double miny,maxy;
        int topScanLine,bottomScanLine;

        miny=max(bottomY,getYMinPoint(t));
        maxy=min(topY,getYMaxPoint(t));

        topScanLine=ceil((topY-maxy)/dy), bottomScanLine=floor((topY-miny)/dy);

        for(int i=topScanLine;i<=bottomScanLine;i++)
        {
            double xa,xb,za,zb,y;
            y=topY- i*dy;

            if(EQ(getYMinPoint(t),getYMaxPoint(t))){
                xa=getXMinPoint(t);
                xb=getXMaxPoint(t);

                za=getZMinPoint(t);
                zb=getZMaxPoint(t);
            }
            else{
                point point0=t.points[0],point1=t.points[1],point2=t.points[2];
                double minx,maxx;
                int leftColInt,rightColInt;
                if(min(point0.y,point1.y)<y && y<max(point0.y,point1.y)){
                    xa=point0.x +(point1.x-point0.x)*(y-point0.y)/(point1.y-point0.y);
                    za=point0.z +(point1.z-point0.z)*(y-point0.y)/(point1.y-point0.y);
                }
                else if(min(point1.y,point2.y)<y && y<max(point1.y,point2.y)){
                    xa=point1.x +(point2.x-point1.x)*(y-point1.y)/(point2.y-point1.y);
                    za=point1.z +(point2.z-point1.z)*(y-point1.y)/(point2.y-point1.y);
                }
                else{
                    xa=point0.x +(point2.x-point0.x)*(y-point0.y)/(point2.y-point0.y);
                    za=point0.z +(point2.z-point0.z)*(y-point0.y)/(point2.y-point0.y);
                }

                if(min(point0.y,point2.y)<y && y<max(point0.y,point2.y)){
                    xb=point0.x +(point2.x-point0.x)*(y-point0.y)/(point2.y-point0.y);
                    zb=point0.z +(point2.z-point0.z)*(y-point0.y)/(point2.y-point0.y);
                }
                else if(min(point1.y,point2.y)<y && y<max(point1.y,point2.y)){
                    xb=point1.x +(point2.x-point1.x)*(y-point1.y)/(point2.y-point1.y);
                    zb=point1.z +(point2.z-point1.z)*(y-point1.y)/(point2.y-point1.y);
                }
                else{
                    xb=point0.x +(point1.x-point0.x)*(y-point0.y)/(point1.y-point0.y);
                    zb=point0.z +(point1.z-point0.z)*(y-point0.y)/(point1.y-point0.y);
                }

                if(xa>xb) {
                    swap(xa,xb);
                    swap(za,zb);
                }

                maxx=min(rightX,xb);
                minx=max(leftX,xa);

                leftColInt=round((minx-leftX)/dx);
                rightColInt=round((maxx-leftX)/dx);

                if(EQ(xa,xb) &&leftColInt==rightColInt){
                    if(za>zb)swap(za,zb);

                    int j=leftColInt;

                    if(((za<zMin) || EQ(za,zMin)) && zMin<zb ){
                        zbuffer[i][j]=zMin +2.0 *__DBL_EPSILON__;
                        image.set_pixel(j,i,t.rgb[0],t.rgb[1],t.rgb[2]);

                    }
                    else if(zMin<za &&za<zbuffer[i][j] && zMin<zb){
                        zbuffer[i][j]=za;
                        //cout<<za<<" ";
                        image.set_pixel(j,i,t.rgb[0],t.rgb[1],t.rgb[2]);
  
                    }
                    continue;
                }
                if(EQ(xa,xb) && leftColInt<rightColInt) continue;
                for(int j=leftColInt;j<=rightColInt;j++){
                    double x,z;
                    x=leftX+j*dx;
                    z=za + (zb-za)*(x-xa)/(xb-xa);
                    if(zMin<z && z<zbuffer[i][j]){
                        zbuffer[i][j]=z;
                        //cout<<z<<" ";
                        image.set_pixel(j,i,t.rgb[0],t.rgb[1],t.rgb[2]);
                    }
                }
            }
        }
    }
    
    // Save the image after processing all triangles
    image.save_image("out.bmp");
    
    // Write z-buffer to file
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            if(zbuffer[i][j]<zMax) {
                o4<<fixed << setprecision(6) << zbuffer[i][j] <<"\t";
                //cout<<zbuffer[i][j]<<" ";
            }
        }
        o4<<"\n";
    }
    
    // Clean up memory
    for(int i=0;i<height;i++) delete[] zbuffer[i];
    delete[] zbuffer;
    
    in.close();
    o4.close();
    
    return 0;
}